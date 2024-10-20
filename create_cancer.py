# This is the code of cancer production
# A code that reproduces repeatedly and at a high speed
#This cancer has the greatest effect on the brain by changing the gene mutation
# And it can also escape from the body's immune system by multiplying a lot
import numpy as np
import random
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def generate_random_dna_sequence(length, mutation_rate):
    nucleotides = ['A', 'C', 'G', 'T']
    sequence = ''.join(random.choice(nucleotides) for _ in range(length))
    mutations = int(length * mutation_rate)
    for _ in range(mutations):
        index = random.randint(0, length - 1)
        sequence = sequence[:index] + random.choice(nucleotides) + sequence[index + 1:]
    return sequence

def replicate_dna(dna_sequence, mutation_rate):
    new_sequence = ''
    for nucleotide in dna_sequence:
        if random.random() < mutation_rate:
            new_sequence += random.choice(['A', 'C', 'G', 'T'])
        else:
            new_sequence += nucleotide
    return new_sequence

def induce_cancerous_mutations(dna_sequence, cancer_mutation_rate, immune_evasion_rate, brain_cancer_factor=1.5):
    new_sequence = ''
    i = 0
    while i < len(dna_sequence):
        # Increase mutation chance for brain-affecting genes
        if random.random() < cancer_mutation_rate * brain_cancer_factor:
            mutation_type = random.choice(['substitution', 'deletion', 'insertion', 'duplication'])
            if mutation_type == 'substitution':
                new_sequence += random.choice(['A', 'C', 'G', 'T'])
                i += 1
            elif mutation_type == 'deletion':
                i += 1  # Skip this nucleotide
            elif mutation_type == 'insertion':
                new_sequence += dna_sequence[i]  # Keep current nucleotide
                new_sequence += random.choice(['A', 'C', 'G', 'T'])  # Insert a random nucleotide
                i += 1
            elif mutation_type == 'duplication':
                new_sequence += dna_sequence[i]  # Keep current nucleotide
                new_sequence += dna_sequence[i]  # Duplicate it
                i += 1
        elif random.random() < immune_evasion_rate:
            # Introduce immune evasion mutations
            new_sequence += random.choice(['A', 'C', 'G', 'T'])  # Randomly mutate a nucleotide
            i += 1
        else:
            new_sequence += dna_sequence[i]
            i += 1
    return new_sequence

def create_chromosome(num_genes, gene_length, mutation_rate):
    chromosome = []
    for _ in range(num_genes):
        gene = generate_random_dna_sequence(gene_length, mutation_rate)
        chromosome.append(gene)
    return chromosome

def replicate_chromosome(chromosome, mutation_rate, cancer_mutation_rate, immune_evasion_rate):
    new_chromosome = []
    for gene in chromosome:
        new_gene = replicate_dna(gene, mutation_rate)
        new_gene = induce_cancerous_mutations(new_gene, cancer_mutation_rate, immune_evasion_rate)
        new_chromosome.append(new_gene)
    return new_chromosome

def reproduce(chromosome, mutation_rate, cancer_mutation_rate, immune_evasion_rate, num_generations, brain_cancer_factor=1.5):
    mutation_data = []  # To store mutation data for visualization
    for generation in range(num_generations):
        # Count mutations in the current generation
        total_mutations = sum(len(gene) - gene.count('A') for gene in chromosome)  # Count non-A nucleotides as mutations
        mutation_data.append({'Generation': generation, 'Total Mutations': total_mutations})

        # Optionally increase mutation rates in certain generations
        if generation % 5 == 0:  # For example, every 5 generations
            mutation_rate *= 1.5  # Increase normal mutation rate
            cancer_mutation_rate *= 6.5  # Increase cancer mutation rate
            immune_evasion_rate *= 4  # Increase immune evasion rate
        
        # Use the brain cancer factor in inducing mutations
        chromosome = replicate_chromosome(chromosome, mutation_rate, cancer_mutation_rate, immune_evasion_rate)
    
    return chromosome, mutation_data

# Parameters
num_genes = 10
gene_length = 100
mutation_rate = 0.01  # Normal mutation rate
cancer_mutation_rate = 1.0  # Initial cancer mutation rate
immune_evasion_rate = 0.9  # Initial immune evasion rate
num_generations = 10

    # Simulation
chromosome = create_chromosome(num_genes, gene_length, mutation_rate)
chromosome, mutation_data = reproduce(chromosome, mutation_rate, cancer_mutation_rate, immune_evasion_rate, num_generations)

    # Convert mutation data to a DataFrame for visualization
mutation_df = pd.DataFrame(mutation_data)

    # Visualization using Seaborn
plt.figure(figsize=(10, 6))
sns.lineplot(data=mutation_df, x='Generation', y='Total Mutations', marker='o')
plt.title('Total Mutations Across Generations')
plt.xlabel('Generation')
plt.ylabel('Total Mutations')
plt.xticks(mutation_df['Generation'])  # Show all generation ticks
plt.grid()
plt.show()
print(chromosome)