import sys
import os
import numpy as np
import pandas as pd

import matplotlib as mpl
mpl.use('QtAgg')  # Use the correct backend for PyQt6
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from PyQt6 import uic
from PyQt6.QtWidgets import (
    QMainWindow, QApplication,
    QLabel, QToolBar, QStatusBar, QCheckBox, QPushButton, QDialog, QDialogButtonBox, QVBoxLayout, QFileDialog,
    QGridLayout, QWidget, QMenu, QHBoxLayout, QTableView, QSizePolicy, QSplitter
)
from PyQt6.QtGui import QAction, QIcon, QKeySequence, QStandardItemModel, QStandardItem
from PyQt6.QtCore import Qt, QAbstractTableModel, QModelIndex

from Model.microbiome_class import MicrobiomeDataAnalyzer
from Model.modificator import otu_table, tax_table

from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=8, height=6, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super().__init__(fig)

class TableModel(QAbstractTableModel):
    def __init__(self, data):
        super().__init__()
        self._data = data

    def rowCount(self, index):
        return self._data.shape[0]

    def columnCount(self, parnet=None):
        return self._data.shape[1]

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if index.isValid():
            if role == Qt.ItemDataRole.DisplayRole or role == Qt.ItemDataRole.EditRole:
                value = self._data.iloc[index.row(), index.column()]
                return str(value)

    def headerData(self, col, orientation, role):
        if orientation == Qt.Orientation.Horizontal and role == Qt.ItemDataRole.DisplayRole:
            return self._data.columns[col]
    def flags(self, index):
        return Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled

class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()
        uic.loadUi("statapp.ui", self)

        # Connect the existing action from the UI file
        self.menuImport.triggered.connect(self.open_file_dialog)
        self.plot_top_button.clicked.connect(self.plot_stat)
        self.plot_rank_button.clicked.connect(self.plot_rank)
        self.plot_pcoa_button.clicked.connect(self.plot_pcoa)
        self.anova_button.clicked.connect(self.anova)
        self.ttest_button.clicked.connect(self.t_test)
        self.textEdit.setText("Welcome to the StatApp! \nPlease import your file for statistical analysis.")

        # setup default plot settings
        mpl.rcParams['font.size'] = 8
        mpl.rcParams['axes.titlesize'] = 8
        mpl.rcParams['axes.labelsize'] = 8
        mpl.rcParams['xtick.labelsize'] = 8
        mpl.rcParams['ytick.labelsize'] = 8
        mpl.rcParams['legend.fontsize'] = 8

# Random forest
# how to disconnect gui from the calculation process
# make on background process

    def open_file_dialog(self):
        dialog = QFileDialog(self)
        dialog.setDirectory(r'C:\Users\egoro\Downloads')
        dialog.setFileMode(QFileDialog.FileMode.ExistingFile)
        dialog.setNameFilter("Tab-sep. files (*.csv *.tsv)")
        dialog.setViewMode(QFileDialog.ViewMode.List)
        if dialog.exec():
            selected_path = dialog.selectedFiles()
            if selected_path:
                filepath = selected_path[0]
                filename = os.path.basename(filepath)
                otu_mat, otu_mat_copy, tax_mat, metadata = self.process_file(filepath)
                self.data_input = MicrobiomeDataAnalyzer(otu_mat, tax_mat, metadata) # Does saving input data as attribute make sense?
                # beta_div = self.data_input.beta_diversity()
        DT = otu_mat_copy
        model = TableModel(DT)
        self.tableView.setModel(model)
        #self.tableView.setVerticalHeaderLabels(all_taxa.to_list())
        self.textEdit.setText(f"File name: {filename}")
        # Connect click signal
        self.tableView.clicked.connect(self.on_taxa_cell_clicked)

    def process_file(self, file_path):
        otu_mat = otu_table(file_path)
        tax_mat = tax_table(file_path)
        last_assigned_taxon = tax_mat.apply(lambda row: row.dropna().iloc[-1] if not row.dropna().empty else np.nan,
                                            axis=1)
        otu_mat_copy = otu_mat.copy()
        otu_mat_copy.insert(0, 'Lowest Taxa', last_assigned_taxon)
        metadata = {
            'SampleID': ['Alice00-1mio.daa', 'Alice01-1mio.daa',
                         'Alice03-1mio.daa', 'Alice06-1mio.daa',
                         'Alice08-1mio.daa', 'Alice34-1mio.daa',
                         'Bob00-1mio.daa', 'Bob01-1mio.daa',
                         'Bob03-1mio.daa', 'Bob06-1mio.daa',
                         'Bob08-1mio.daa', 'Bob34-1mio.daa'],
            'Group': ['Alice', 'Alice', 'Alice',
                      'Alice', 'Alice', 'Alice',
                      'Bob', 'Bob', 'Bob',
                      'Bob', 'Bob', 'Bob'],
            'Property': ['0-', '1+', '3+', '6+', '8-', '34-',
                         '0-', '1+', '3+', '6+', '8-', '34-']
        }
        metadata_alice_bob = pd.DataFrame(metadata)
        return otu_mat, otu_mat_copy, tax_mat, metadata_alice_bob

    # Cleans subplots for new plotting functions
    def clear_canvas(self):
        self.graph_widget.axes.clear()

    def t_test(self):
        self.textEdit.clear()
        result = self.data_input.t_test()
        # Display the results
        for otu, (t_stat, p_value) in result.items():
            self.textEdit.append(f"{otu}: T statistic = {t_stat}, P-value = {p_value}")

    def anova(self):
        self.textEdit.clear()
        result = self.data_input.anova_test()
        # Display the results
        for otu, (f_value, p_value) in result.items():
            self.textEdit.append(f"{otu}: F-value = {f_value}, P-value = {p_value}")

    def plot_stat(self):
        self.clear_canvas()
        # Use the input's plot_top function to render the plot
        top_plot, top_df = self.data_input.plot_top(5, self.graph_widget)

        # Connect the event handler to the canvas
        self.graph_widget.mpl_connect("button_press_event", lambda event: self.on_click(event, top_plot, top_df))

    def plot_rank(self):
        self.clear_canvas()

        # Use the input's plot_top function to render the plot
        rank_plot, rank_df = self.data_input.plot_rank("Phylum", self.graph_widget)

        # Connect the event handler to the canvas
        self.graph_widget.mpl_connect("button_press_event", lambda event: self.on_click(event, rank_plot, rank_df))

    def plot_pcoa(self):
        self.clear_canvas()
        self.data_input.plot_pcoa(self.graph_widget)

    def on_click(self, event, plot, df):
        totals = df.sum(axis=1)
        # Check if the click occurred inside the axes
        if event.inaxes == self.graph_widget.axes:
            for bar_group, group_name in zip(plot.containers, df.columns):
                for rect, value, total in zip(bar_group, df[group_name], totals):
                    height = rect.get_height()
                    # Check if the mouse click is inside the rectangle
                    if rect.contains(event)[0]:  # Returns (True, details) if inside
                        taxon = group_name
                        percentage = f'{(value / total) * 100:.1f}%'  # Calculate percentage
                        self.statusBar().showMessage(f"Taxon:{taxon}, Percentage: {percentage}, Read Count: {height}")
                        break

    def on_plot_click(self, event):
        if event.inaxes:  # Check if the click is inside the plot area
            x, y = event.xdata, event.ydata
            # Display the clicked coordinates
            self.statusBar().showMessage(f"Clicked at: x={x:.2f}, y={y:.2f}")

    def on_taxa_cell_clicked(self, index: QModelIndex):
        from ete3 import NCBITaxa
        ncbi = NCBITaxa()
        # index holds both row & column
        if index.column() != 0:
            return

        # Retrieve the taxon name
        taxa_name = index.data()

        # Lookup full classification
        try:
            tax_tab = self.data_input.Taxa_table
            # Find location of the value
            # matches = tax_tab[tax_tab == taxa_name]
            # Step 1: Convert name to taxid
            name2taxid = ncbi.get_name_translator([taxa_name])
            taxid = name2taxid[taxa_name][0]

            # Step 2: Get full lineage taxids
            lineage = ncbi.get_lineage(taxid)

            # Step 4: (Optional) Get ranks for each taxid
            major_ranks = ["superkingdom", "domain", "phylum", "class", "order", "family", "genus", "species"]
            ranks = ncbi.get_rank(lineage)
            sorted_ranks = {k: v for (k, v) in ranks.items() if v in major_ranks}
            sorted_lineage = [k for (k,v) in sorted_ranks.items()]
            # Step 3: Map lineage taxids to names
            names = ncbi.get_taxid_translator(sorted_lineage)
            # Build ordered path (Domain → ... → species/order/etc.)

            taxonomy_path = [names[t] for t in sorted_lineage]

        except KeyError:
            self.statusBar().showMessage(f"No full taxonomy found for '{taxa_name}'")
            return

        # Format for display

        self.statusBar().showMessage(', '.join(map(str, taxonomy_path)))

if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = MainWindow()
    w.show()
    sys.exit(app.exec())

# ToDo:
# 1. Add descriptions to OTUs in the table (Drop down).
# Or: lowest rank in the table, once user clicks on the cell, show whole phylogeny on the left
# 2. Give a choice: relative or absolute abundancies
# 3. Create a separate table with test results?
# 4. For testing: for each row in the newly created table show a plot
# 5. Plot rank: add other rank possibilities
# Keep everything basic, look for traditional GUI interface layout.