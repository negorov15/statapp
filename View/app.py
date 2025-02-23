import sys
import datatable as dt
import numpy as np
import pandas as pd
from pathlib import Path
import plotly.express as px

import matplotlib as mpl
mpl.use('QtAgg')  # Use the correct backend for PyQt6
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg,
    NavigationToolbar2QT as NavigationToolbar,
)

from PyQt6 import uic
from PyQt6.QtWidgets import (
    QMainWindow, QApplication,
    QLabel, QToolBar, QStatusBar, QCheckBox, QPushButton, QDialog, QDialogButtonBox, QVBoxLayout, QFileDialog,
    QGridLayout, QWidget, QMenu, QHBoxLayout, QTableView, QSizePolicy, QSplitter
)
from PyQt6.QtGui import QAction, QIcon, QKeySequence, QStandardItemModel, QStandardItem
from PyQt6.QtCore import Qt, QSize, QAbstractTableModel

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
            file_name = dialog.selectedFiles()
            if file_name:
                file_path = file_name[0]
                otu_mat, tax_mat, metadata = self.process_file(file_path)
                self.data_input = MicrobiomeDataAnalyzer(otu_mat, tax_mat, metadata) # Does saving input data as attribute make sense?
        DT = self.data_input.OTU_table
        model = TableModel(DT)
        self.tableView.setModel(model)

    def process_file(self, file_path):
        otu_mat = otu_table(file_path)
        tax_mat = tax_table(file_path)
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
        return otu_mat, tax_mat, metadata_alice_bob

    # Cleans subplots for new plotting functions
    def clear_canvas(self):
        self.graph_widget.axes.clear()

    def plot_stat(self):
        self.clear_canvas()
        # Use the input's plot_top function to render the plot
        top_plot, top_df = self.data_input.plot_top(5, self.graph_widget)
        # Calculate row-wise totals for percentages
        totals = top_df.sum(axis=1)

        def on_click(event):
            # Check if the click occurred inside the axes
            if event.inaxes == self.graph_widget.axes:
                for bar_group, group_name in zip(top_plot.containers, top_df.columns):
                    for rect, value, total in zip(bar_group, top_df[group_name], totals):
                        height = rect.get_height()
                        # Check if the mouse click is inside the rectangle
                        if rect.contains(event)[0]:  # Returns (True, details) if inside
                            taxon = group_name
                            percentage = f'{(value / total) * 100:.1f}%'  # Calculate percentage
                            self.statusBar().showMessage(f"Taxon:{taxon}, Percentage: {percentage}, Read Count: {height}")
                            break

        # Connect the event handler to the canvas
        self.graph_widget.mpl_connect("button_press_event", on_click)

    def plot_rank(self):
        self.clear_canvas()

        # Use the input's plot_top function to render the plot
        rank_plot, rank_df = self.data_input.plot_rank("Phylum", self.graph_widget)
        # Calculate row-wise totals for percentages
        totals = rank_df.sum(axis=1)
        def on_click(event):
            # Check if the click occurred inside the axes
            if event.inaxes == self.graph_widget.axes:
                for bar_group, group_name in zip(rank_plot.containers, rank_df.columns):
                    for rect, value, total in zip(bar_group, rank_df[group_name], totals):
                        height = rect.get_height()
                        # Check if the mouse click is inside the rectangle
                        if rect.contains(event)[0]:  # Returns (True, details) if inside
                            taxon = group_name
                            percentage = f'{(value / total) * 100:.1f}%'  # Calculate percentage
                            self.statusBar().showMessage(f"Taxon:{taxon}, Percentage: {percentage}, Read Count: {height}")
                            break

        # Connect the event handler to the canvas
        self.graph_widget.mpl_connect("button_press_event", on_click)

    def plot_pcoa(self):
        ...

    def on_plot_click(self, event):
        if event.inaxes:  # Check if the click is inside the plot area
            x, y = event.xdata, event.ydata
            # Display the clicked coordinates
            self.statusBar().showMessage(f"Clicked at: x={x:.2f}, y={y:.2f}")

app = QApplication(sys.argv)
w = MainWindow()
w.show()
app.exec()