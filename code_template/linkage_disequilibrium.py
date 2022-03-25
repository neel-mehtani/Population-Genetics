import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def split_phased_snps(phased_snps):
    """Split phased SNP data into sequential pair of maternal and paternal SNPs.

    Arguments:
        phased_snps (pd.DataFrame): a (num_snps x (9 + num_samples)) dataframe loaded from a VCF file, with first
            nine columns ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] and the
            following columns containing the sequenced genotypes, where each value is a string of the form

                                    "maternal_allele | paternal_allele"

            for example, an entry "0 | 1" indicates that the maternal chromosome contained the reference allele
            at that locus while the paternal chromosome contained the alternate allele.

    Returns:
        split_snps: a (2 x num_snps x num_samples)-dimension Numpy array of strings "0" or "1"s, where snps[0] is the
            maternal haplotype for all samples and snps[1] is the paternal haplotype for all samples
    """
    pass

def compute_allele_frequencies(snps):
    """Compute reference and alternate allele frequencies for each SNP.

    Arguments:
        snps: a (num_snps x (2 * num_samples))-dimension Numpy array of strings "0" or "1"s, where the first num_samples
            columns are alleles from the maternal chromosome and the second num_samples columns are the corresponding
            alleles from the paternal chromosome.

    Returns:
        a 2-tuple of (snp_ref_frequency, snp_alt_frequency) representing the reference and alternate frequencies
            for each snp respectively, where snp_ref_frequency[i] + snp_alt_frequency[i] = 1
    """
    pass

def compute_haplotype_frequencies(snps):
    """Compute haplotype frequencies for every pair of SNPs.

    Arguments:
        snps: a (num_snps x (2 * num_samples))-dimension Numpy array of strings "0" or "1"s, where the first num_samples
            columns are alleles from the maternal chromosome and the second num_samples columns are the corresponding
            alleles from the paternal chromosome.

    Returns:
        a 4-tuple of (p_11, p_12, p_21, p_22) representing the haplotype frequencies of ("0", "0"), ("0", "1"),
            ("1", "0") and ("1", "1") for each pair respectively, where each element is a (num_snps x num_snps) Numpy
            array, and where the sum of p_11, p_12, p_21, p_22 for each pair of snps is 1.
    """
    pass

def calculate_D(haplotype_frequencies):
    """Calculate the linkage disequilibrium D.

    Arguments:
        haplotype_frequencies: a 4-tuple of (p_11, p_12, p_21, p_22) representing the haplotype frequencies of ("0", "0"), ("0", "1"),
            ("1", "0") and ("1", "1") for each pair respectively, where each element is a (num_snps x num_snps) Numpy
            array, and where the sum of p_11, p_12, p_21, p_22 for each pair of snps is 1.

    Returns:
        a (num_snps x num_snps) Numpy array of the raw linkage disequilibrium estimates for each pair of snps. This ranges in value
            from -0.25 to 0.25.
    """
    pass

def calculate_D_prime(allele_frequencies, haplotype_frequencies):
    """Calculate the standardized linkage disquilibrium D'.

    Arguments:
        allele_frequencies: a 2-tuple of (snp_ref_frequency, snp_alt_frequency) representing the reference and alternate frequencies
            for each allele respectively, where snp_ref_frequency[i] + snp_alt_frequency[i] = 1
        haplotype_frequencies: a 4-tuple of (p_11, p_12, p_21, p_22) representing the haplotype frequencies of ("0", "0"), ("0", "1"),
            ("1", "0") and ("1", "1") for each pair respectively, where each element is a (num_snps x num_snps) Numpy
            array, and where the sum of p_11, p_12, p_21, p_22 for each pair of snps is 1.

    Returns:
        a (num_snps x num_snps) Numpy array of the standardized linkage disequilibrium estimates for each pair of snps. This
            ranges in value from 0 to 1.
    """
    pass

def calculate_r_squared(allele_frequencies, haplotype_frequencies):
    """Calculate the square of Pearson's correlation coefficient r^2.

    Arguments:
        allele_frequencies: a 2-tuple of (snp_ref_frequency, snp_alt_frequency) representing the reference and alternate frequencies
            for each allele respectively, where snp_ref_frequency[i] + snp_alt_frequency[i] = 1
        haplotype_frequencies: a 4-tuple of (p_11, p_12, p_21, p_22) representing the haplotype frequencies of ("0", "0"), ("0", "1"),
            ("1", "0") and ("1", "1") for each pair respectively, where each element is a (num_snps x num_snps) Numpy
            array, and where the sum of p_11, p_12, p_21, p_22 for each pair of snps is 1.

    Returns:
        a (num_snps x num_snps) Numpy array of the r^2 values for each pair of snps. This
            ranges in value from 0 to 1.
    """
    pass

if __name__ == "__main__":
    phased_snps = pd.read_csv("../provided_data/chromosome_1_phased_first_5000_snps.vcf", delimiter="\t")
    # TODO: make calls to above functions
