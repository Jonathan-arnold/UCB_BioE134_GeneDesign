from genedesign.seq_utils.reverse_complement import reverse_complement


class InternalRBSChecker:
    """
    Checks for internal ribosome binding sites (Shine-Dalgarno-like sequences)
    upstream of ATG codons that could cause unintended translation initiation.
    """

    def __init__(self):
        self.sd_sequences = []

    def initiate(self):
        # Shine-Dalgarno consensus is AAGGAG; include common variants
        # Ordered from strongest/longest to shortest match
        self.sd_sequences = [
            "AAGGAG",  # full consensus
            "AGGAGG",  # extended core
            "AAGGAG",  # full match
            "AGGAG",   # strong core
            "AAGGA",   # 5' partial
            "GGAGG",   # 3' extended
            "AGGA",    # minimal strong
            "GGAG",    # minimal core
            "GAGG",    # minimal variant
        ]
        # Deduplicate while preserving order
        seen = set()
        deduped = []
        for s in self.sd_sequences:
            if s not in seen:
                seen.add(s)
                deduped.append(s)
        self.sd_sequences = deduped

    def run(self, dnaseq):
        """
        Scans for Shine-Dalgarno-like motifs 4-12 bp upstream of any ATG codon.

        Parameters:
            dnaseq (str): DNA sequence to check.

        Returns:
            tuple: (True, None) if no internal RBS found,
                   (False, str) with the offending SD+ATG region if found.
        """
        seq = dnaseq.upper()
        rc = reverse_complement(dnaseq).upper()
        combined = seq + "x" + rc

        # Find all ATG positions in the combined sequence
        for i in range(len(combined) - 2):
            if combined[i:i+3] != "ATG":
                continue

            # Check the region 4-12 bp upstream of this ATG for SD motifs
            # The SD sequence is typically spaced 4-12 nt before the start codon
            upstream_start = max(0, i - 12)
            upstream_end = max(0, i - 4)
            if upstream_end <= upstream_start:
                continue

            upstream = combined[upstream_start:upstream_end]

            for sd in self.sd_sequences:
                if sd in upstream:
                    # Return the full context: SD region + spacing + ATG
                    context_start = max(0, i - 12)
                    context = combined[context_start:i+3]
                    return False, context

        return True, None


def main():
    checker = InternalRBSChecker()
    checker.initiate()
    # Should fail: strong SD sequence upstream of ATG
    result = checker.run("AAAAAAGGAGATCGATG")
    print(result)
    # Should pass: no SD near ATG
    result = checker.run("TTTCCCGGGAAATTTCCC")
    print(result)


if __name__ == "__main__":
    main()
