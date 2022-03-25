import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats.distributions import chi2
from sklearn.preprocessing import StandardScaler

def snp_array(vcf_df, numSnps, numSamples):
    
    split_snps = np.ndarray(shape=(2, numSnps, numSamples))
    snps = vcf_df.iloc[:, 9:]
    
    for i in range(numSnps):
        for j in range(numSamples):
            
            split_snps[0, i, j] = snps.iloc[i, j][0]
            split_snps[1, i, j] = snps.iloc[i, j][-1]
    
    return split_snps

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
        snps: a (2 x num_snps x num_samples)-dimension Numpy array of strings "0" or "1"s, where snps[0] is the
            maternal haplotype for all samples and snps[1] is the paternal haplotype for all samples
    """
    # TODO: copy your implementation of this function from linkage_disequilibrium.py
    
    data = []
    
    with open(vcf_file, 'r') as invcf:

        for line in invcf:
            print(ct)
            if line.startswith('#'):
                continue
            
            line = line.strip().split()
            data.append(line)
            
    
    numSamples = len(data[0]) - 9
    
    vcf_cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    vcf_cols.extend(["obs_{}".format(i+1) for i in range(numSamples)])
    
    res_df = pd.DataFrame(data, columns=vcf_cols)
    
    return res_df
    

def compute_effective_allele_frequencies(split_snps):
    """Compute reference and alternate allele frequencies for each SNP.

    Arguments:
        split_snps: a (2 x num_snps x num_samples)-dimension Numpy array of strings "0" or "1"s, where snps[0] is the
            maternal haplotype for all samples and snps[1] is the paternal haplotype for all samples
    
    Returns:
        a 2-tuple of (snp_ref_frequency, snp_alt_frequency) representing the reference and alternate frequencies
            for each snp respectively, where snp_ref_frequency[i] + snp_alt_frequency[i] = 1
    """
    genotype_cts = compute_genotype_counts(split_snps)
    numSnps = len(split_snps[0])
    ref_freq, alt_freq = [], []     
    
    for i in range(numSnps):
        
        total_cts = genotype_cts[0][i] + genotype_cts[1][i] + genotype_cts[2][i]
        curr_ref_freq = (2 * genotype_cts[0][i] + (genotype_cts[1][i]))/(2*total_cts)
        curr_alt_freq = (2 * genotype_cts[2][i] + (genotype_cts[1][i]))/(2*total_cts)
        
        ref_freq.append(curr_ref_freq)
        alt_freq.append(curr_alt_freq)
        
    allele_freqs = (np.array(ref_freq), np.array(alt_freq))
    
    return allele_freqs

def compute_genotype_counts(split_snps):
    """Compute reference and alternate allele frequencies for each SNP.

    Arguments:
        split_snps: a (2 x num_snps x num_samples)-dimension Numpy array of strings "0" or "1"s, where snps[0] is the
            maternal haplotype for all samples and snps[1] is the paternal haplotype for all samples
    
    Returns:
        a 3-tuple of (homozygous_reference_counts, heterozygous_counts, homozygous_alternate_counts) representing the
            homozygous reference, heterozygous and homozygous alternate counts for each snp respectively,
            where snp_ref_frequency[i] + snp_alt_frequency[i] = 1
    """
    homo_ref_cts, hetero_cts, homo_alt_cts = [], [], []
    numSnps = len(split_snps[0])
    
    for i in range(numSnps):
        curr_genotypes = list(zip(split_snps[0][i], split_snps[1][i]))
        cts = dict(Counter(curr_genotypes))
        homo_ref_cts.append(cts[(0.0, 0.0)])
        homo_alt_cts.append(cts[(1.0, 1.0)])
        hetero_cts.append(cts[(0.0, 1.0)] + cts[(1.0, 0.0)])
    
    genotype_counts = (np.array(homo_ref_cts), np.array(hetero_cts), np.array(homo_alt_cts))
    
    return genotype_counts

def calculate_chi_squared_statistic(effective_allele_frequencies, genotype_counts):
    """Calculate the χ2 statistic for each SNP.

    Arguments:
        effective_allele_frequencies: a 2-tuple of (snp_ref_frequency, snp_alt_frequency) representing the reference and alternate frequencies
            for each snp respectively, where snp_ref_frequency[i] + snp_alt_frequency[i] = 1
        genotype_counts: a 3-tuple of (homozygous_reference_counts, heterozygous_counts, homozygous_alternate_counts)
            representing the homozygous reference, heterozygous and homozygous alternate counts for each snp respectively,
            where snp_ref_frequency[i] + snp_alt_frequency[i] = 1
    
    Returns:
        a num_snps-dimension Numpy array containing the χ2 statistic for each SNP.
                
    """
    numSnps = len(genotype_counts[0])
    numSamples = genotype_counts[0][0] + genotype_counts[1][0] + genotype_counts[2][0]
    
    expected_genotype_counts = [np.zeros(numSnps), np.zeros(numSnps), np.zeros(numSnps)]
    
    for i in range(numSnps):
        
        expected_genotype_counts[0][i] = numSamples*(effective_allele_frequencies[0][i]**2)
        expected_genotype_counts[1][i] = numSamples*(2*effective_allele_frequencies[0][i]*effective_allele_frequencies[1][i])
        expected_genotype_counts[2][i] = numSamples*(effective_allele_frequencies[1][i]**2)
    
    chi_stats = ((np.array(expected_genotype_counts) - np.array(genotype_counts))**2)/np.array(expected_genotype_counts)
    chi_stats = np.sum(chi_stats, axis = 0)
            
    return chi_stats

def detect_snps_under_selection(chi_squared_statistic, alpha=1e-2, dof=1):
    """Determie which SNPs violate Hardy-Weinberg equilibrium according to the χ2 statistic.

    Arguments:
        a num_snps-dimension Numpy array containing the χ2 statistic for each SNP.
    
    Returns:
        a num_snps-dimension boolean Numpy array representing whether each SNP violates or doesn't violate HW-equilibrium
    """

    is_eq = []
    p_vals = []
    
    for i in chi_squared_statistic:
        
        p_value = 1 - chi2.cdf(i, dof)
        p_vals.append(p_value)
        
        if p_value < alpha:
            is_eq.append(True)
        else:
            is_eq.append(False)

    return np.array(is_eq, dtype=bool), p_vals

if __name__ == "__main__":
    phased_snps = pd.read_csv("../provided_data/chromosome_19_phased_snps.vcf", delimiter="\t")
    
