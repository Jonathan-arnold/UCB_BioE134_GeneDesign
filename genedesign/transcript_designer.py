import random
import os
from collections import defaultdict
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.codon_checker import CodonChecker


class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence and chooses an RBS
    using weighted random codon selection based on E. coli codon usage frequencies.
    """

    def __init__(self):
        self.codonUsage = {}  # aa -> ([codons], [weights])
        self.rbsChooser = None
        self.forbiddenChecker = None
        self.promoterChecker = None
        self.codonChecker = None

    def initiate(self) -> None:
        """
        Initializes the codon usage table from codon_usage.txt and the RBS chooser.
        """
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        self.forbiddenChecker = ForbiddenSequenceChecker()
        self.forbiddenChecker.initiate()

        self.promoterChecker = PromoterChecker()
        self.promoterChecker.initiate()

        self.codonChecker = CodonChecker()
        self.codonChecker.initiate()

        # Parse codon usage file
        usage_file = os.path.join(os.path.dirname(__file__), "data", "codon_usage.txt")
        aa_to_codons = defaultdict(list)
        aa_to_weights = defaultdict(list)

        with open(usage_file, "r") as f:
            for line in f:
                parts = line.split()
                if len(parts) < 3:
                    continue
                codon = parts[0]
                aa = parts[1]
                fraction = float(parts[2])
                if aa == "*":  # skip stop codons
                    continue
                aa_to_codons[aa].append(codon)
                aa_to_weights[aa].append(fraction)

        self.codonUsage = {
            aa: (aa_to_codons[aa], aa_to_weights[aa]) for aa in aa_to_codons
        }

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates the peptide sequence to DNA and selects an RBS.

        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.

        Returns:
            Transcript: The transcript object with the selected RBS and translated codons.
        """
        max_attempts = 1000
        window_size = 3  # number of new codons to lock in per window
        lookahead = 6  # number of lookahead codons

        selectedRBS = self.rbsChooser.run("", ignores)
        utr = selectedRBS.utr.upper()

        codons = []
        total_rare_codons = 0

        i = 0  # current position in peptide
        while i < len(peptide):
            remaining = len(peptide) - i
            n_select = min(window_size, remaining)
            n_lookahead = min(lookahead, remaining - n_select)

            found = False

            for _ in range(max_attempts):
                candidate = []
                for j in range(n_select + n_lookahead):
                    aa = peptide[i + j]
                    options, weights = self.codonUsage[aa]
                    candidate.append(random.choices(options, weights=weights)[0])

                # Build the window DNA: prior context + candidate
                # Use enough context so window_dna >= 29bp (promoter checker's sliding frame)
                if i == 0:
                    prior_dna = utr
                else:
                    prior_dna = "".join(codons[max(0, len(codons) - 10) :])

                candidate_dna = "".join(candidate)
                if i + n_select >= len(peptide):
                    candidate_dna += "TAA"

                window_dna = prior_dna + candidate_dna

                # Check 1: forbidden sequences
                passed_forbidden, _ = self.forbiddenChecker.run(window_dna)
                if not passed_forbidden:
                    continue

                # Check 2: internal promoters
                passed_promoter, _ = self.promoterChecker.run(window_dna)
                if not passed_promoter:
                    continue

                # Check 3: codon usage quality on the local window
                codon_window = codons + candidate
                _, diversity, _, cai = self.codonChecker.run(codon_window)

                # Count rare codons in the candidate selection only
                candidate_rare = sum(
                    1
                    for c in candidate[:n_select]
                    if c in self.codonChecker.rare_codons
                )

                # Determine allowed rare codon budget based on position
                end_pos = i + n_select
                if end_pos <= len(peptide) / 3:
                    allowed_rare = 1
                elif end_pos <= 2 * len(peptide) / 3:
                    allowed_rare = 2
                else:
                    allowed_rare = 3

                min_diversity = 0.5 * (end_pos / len(peptide))
                if (
                    diversity < min_diversity
                    or (total_rare_codons + candidate_rare) > allowed_rare
                    or cai < 0.2
                ):
                    continue

                # Passed all checks
                total_rare_codons += candidate_rare
                codons.extend(candidate[:n_select])
                found = True
                break

            if not found:
                raise ValueError(
                    f"Failed to find a valid codon window at position {i} "
                    f"after {max_attempts} attempts."
                )

            i += n_select

        # Append the stop codon
        codons.append("TAA")

        # Return the Transcript object
        return Transcript(selectedRBS, peptide, codons)


if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    peptide = "MYPFIRTARMTV"

    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)

    # Print out the transcript information
    print(transcript)
