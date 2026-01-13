import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon
from skbio.stats.distance import permanova
from Model.get_lineage import *
from skbio.stats.ordination import pcoa
from skbio.diversity import beta_diversity
from Model.modificator import otu_table, tax_table

class MicrobiomeDataAnalyzer:

    def __init__(self, OTU_table: str = "", Taxa_table: str = "", Metadata: str = "") -> None:
        self.OTU_table = OTU_table
        self.Taxa_table = Taxa_table
        self.Metadata = Metadata

    def beta_diversity(self):
        dropped_df = self.OTU_table.dropna()
        otumat_transposed = dropped_df.T
        data = otumat_transposed.values
        ids = otumat_transposed.index
        b_div = beta_diversity(metric="braycurtis", counts=data, ids=ids)
        return b_div

    def plot_rank(self, rank, canvas):
        taxa_mapping = self.Taxa_table[rank]
        rank_df = self.OTU_table.groupby(taxa_mapping, axis=0).sum().T

        resulting_plot = rank_df.plot(kind='bar', stacked=True, ax=canvas.axes)

        canvas.axes.set_title(f'Abundance of {rank}')
        canvas.axes.set_ylabel('Abundance')
        canvas.axes.tick_params(axis='x', rotation=45)  # Rotate labels
        canvas.axes.legend(title=rank, loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
        canvas.figure.tight_layout()  # Automatically adjust layout to fit labels
        canvas.draw()

        return resulting_plot, rank_df

    def plot_top(self, top, canvas):
        # Calculate top taxa and plot
        total_counts = self.OTU_table.sum(axis=1)
        top_otus = total_counts.nlargest(top).index

        top_otus_tax_rank = {}
        for otu in top_otus:
            tax_info = self.Taxa_table.loc[otu]
            specific_rank = tax_info.dropna().iloc[-1] if not tax_info.dropna().empty else 'Unknown'
            top_otus_tax_rank[otu] = specific_rank

        filtered_df = self.OTU_table.loc[top_otus]
        filtered_df.index = [top_otus_tax_rank[otu] for otu in filtered_df.index]
        df_t = filtered_df.T
        totals = df_t.sum(axis=1)

        # Plot stacked bar chart
        # pandas' built-in plotting function, which wraps around Matplotlib.
        resulting_plot = df_t.plot(kind='bar', stacked=True, ax=canvas.axes)

        # Customize axes
        canvas.axes.set_title(f'Top {top} taxa by read count')
        canvas.axes.set_ylabel('Read Count')
        # Update the legend (if required)
        canvas.axes.legend(title="Taxa", loc='upper left')

        # Adjust x-axis labels
        canvas.axes.tick_params(axis='x', rotation=45)  # Rotate labels
        canvas.figure.tight_layout()  # Automatically adjust layout to fit labels
        canvas.draw()

        return resulting_plot, df_t

    def plot_pcoa(self, canvas):
        distance = self.beta_diversity()
        pcoa_results = pcoa(distance)
        pcoa_samples = pcoa_results.samples

        groups = np.array(self.Metadata['Group'])
        treatment_day = np.array(self.Metadata['Property'])
        pcoa_samples.insert(0, "Group", groups)
        pcoa_samples.insert(1, "Property", treatment_day)

        unique_groups = pcoa_samples['Group'].unique()
        colors = plt.get_cmap('tab10')(range(len(unique_groups)))
        color_map = dict(zip(unique_groups, colors))

        seen_groups = set()

        for idx, row in pcoa_samples.iterrows():
            group = row['Group']
            pc1 = row['PC1']
            pc2 = row['PC2']
            property = row['Property']

            if group not in seen_groups:
                canvas.axes.scatter(pc1, pc2, color=color_map[group], label=group, alpha=0.7)
                seen_groups.add(group)
            else:
                canvas.axes.scatter(pc1, pc2, color=color_map[group], alpha=0.7)

            canvas.axes.text(pc1, pc2, property, fontsize=8)

        canvas.axes.set_xlabel('PC1', fontsize=14)
        canvas.axes.set_ylabel('PC2', fontsize=14)
        canvas.axes.set_title('PCoA of Beta Diversity', fontsize=16)
        canvas.axes.legend(loc='upper left', bbox_to_anchor=(1, 1))
        canvas.figure.tight_layout()  # Automatically adjust layout to fit labels
        canvas.draw()


    def prepare_dataset(self):
        otumat_transposed = self.OTU_table.T
        otumat_transposed = otumat_transposed.reset_index()
        otumat_transposed = otumat_transposed.rename(columns={'index': 'SampleID'})
        filtered_metadata = self.Metadata.drop("Property", axis='columns')
        input_data = pd.merge(otumat_transposed, filtered_metadata, on='SampleID')
        column_to_move = input_data.pop("Group")
        input_data.insert(1, "Group", column_to_move)
        return input_data

#     # Statistical analysis functions:
#     # 1. t_test()
#     # 2. anova_test()
#     # 3. kruskal()
#     # 4. permanova()
#     # 5. beta_diversity()
#     # Perform an independent t-Test

    def t_test(self):
        # Prepare the dataset
        input_data = self.prepare_dataset()

        # Find all unique color groups in the dataset.
        unique_groups = input_data['Group'].unique()
        # Save the possible colors
        column_property = input_data['Group']

        # Perform t-Test for each OTU column
        t_test_results = {}
        for otu in input_data.columns[2:]:  # Skipping 'SampleID' and 'Property'

            # Filter the dataframe by color and drop missing values.
            groups = [input_data[column_property == property][otu].dropna() for property in unique_groups]

            t_stat, p_value = stats.ttest_ind(*groups)

            # Store the results
            t_test_results[otu] = (t_stat, p_value)

        return t_test_results

    # Perform Anova-Test
    def anova_test(self):
        # Prepare the dataset
        input_data = self.prepare_dataset()
        # Prepare groups for ANOVA
        # Find all unique groups in the dataset.
        unique_groups = input_data['Group'].unique()
        # Save the possible properties (for example: fur color of the subject)
        column_property = input_data['Group']

        # Perform ANOVA for each OTU column
        anova_results = {}
        for otu in input_data.columns[2:]:  # Skipping 'SampleID' and 'Property'

            # Filter the dataframe by color and drop missing values.
            groups = [input_data[column_property == property][otu].dropna() for property in unique_groups]

            f_value, p_value = stats.f_oneway(*groups)

            # Store the results
            anova_results[otu] = (f_value, p_value)

        return anova_results

    # Perform Kruskal test. Alternative to the ANOVA-test. Data does not have to be normally distributed.
    def kruskal(self):
        # Prepare the dataset
        input_data = self.prepare_dataset()
        # Find all unique groups in the dataset.
        unique_groups = input_data['Group'].unique()
        # Save the possible properties (for example: fur color of the subject)
        column_property = input_data['Group']

        # Perform Kruskal for each OTU column
        kruskal_results = {}
        for otu in input_data.columns[2:]:  # Skipping 'SampleID' and 'Property'

            # Filter the dataframe by color and drop missing values.
            groups = [input_data[column_property == property][otu].dropna() for property in unique_groups]

            stat, p_value = stats.kruskal(*groups)

            # Store the results
            kruskal_results[otu] = (stat, p_value)

        # Display the results
        for otu, (stat, p_value) in kruskal_results.items():
            print(f"{otu}: H-statistic = {stat}, P-value = {p_value}")

    def wilcoxon_test(self):
        # Prepare the dataset
        input_data = self.prepare_dataset()

        # Assuming two unique groups represent paired data (e.g., "before" and "after")
        if len(input_data['Group'].unique()) != 2:
            raise ValueError("There must be exactly two groups for the Wilcoxon signed-rank test.")

        # Retrieve the unique groups for clarity
        group1, group2 = input_data['Group'].unique()

        # Perform Wilcoxon test for each OTU column
        wilcoxon_results = {}

        for otu in input_data.columns[2:]:  # Skipping 'SampleID' and 'Group'

            # Filter the dataframe by group and drop missing values
            group1_data = input_data[input_data['Group'] == group1][otu].dropna()
            group2_data = input_data[input_data['Group'] == group2][otu].dropna()

            # Ensure equal lengths before performing the test
            min_len = min(len(group1_data), len(group2_data))
            if min_len > 0:
                stat, p_value = wilcoxon(group1_data[:min_len], group2_data[:min_len])

                # Store the results
                wilcoxon_results[otu] = (stat, p_value)
            else:
                wilcoxon_results[otu] = (None, None)

        # Display the results
        for otu, (stat, p_value) in wilcoxon_results.items():
            print(f"{otu}: Wilcoxon statistic = {stat}, P-value = {p_value}")

    def permanova(self):
        distance = self.beta_diversity()
        grouping = self.Metadata['Group']
        permanova_results = permanova(distance, grouping)
        return permanova_results