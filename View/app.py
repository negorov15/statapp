import sys
import datatable as dt
import numpy as np
import pandas as pd
from pathlib import Path

import matplotlib
matplotlib.use('QtAgg')  # Use the correct backend for PyQt6
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas, FigureCanvasQTAgg

from PyQt6.QtWidgets import (
    QMainWindow, QApplication,
    QLabel, QToolBar, QStatusBar, QCheckBox, QPushButton, QDialog, QDialogButtonBox, QVBoxLayout, QFileDialog,
    QGridLayout, QWidget, QMenu, QHBoxLayout, QTableView, QSizePolicy
)
from PyQt6.QtGui import QAction, QIcon, QKeySequence, QStandardItemModel, QStandardItem
from PyQt6.QtCore import Qt, QSize

from Model.microbiome_class import MicrobiomeDataAnalyzer
from Model.modificator import otu_table, tax_table


class CustomDialog(QDialog):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Import your file!")

        QBtn = QDialogButtonBox.StandardButton.Ok

        self.buttonBox = QDialogButtonBox(QBtn)
        self.buttonBox.accepted.connect(self.accept)

        self.layout = QVBoxLayout()
        message = QLabel("Implement the file open dialog box!")
        self.layout.addWidget(message)
        self.layout.addWidget(self.buttonBox)
        self.setLayout(self.layout)

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=40, height=24, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)

class MainWindow(QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()

        self.setWindowTitle("Microbiome Data Analyzer")
        self.setMinimumSize(QSize(600, 400))

        self.file_path = None
        self.label = None
        self.table_view = None
        self.setup_toolbar()
        self.setup_menu()
        self.setup_layout()

    # Main layout of the application
    def setup_layout(self):
        self.label = QLabel("Statistical Analysis of data!", self)
        self.label.setAlignment(Qt.AlignmentFlag.AlignCenter)

        layout = QVBoxLayout()
        layout.addWidget(self.label)

        central_widget = QWidget()
        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)

    def setup_toolbar(self):
        toolbar = QToolBar("All instruments")
        toolbar.setIconSize(QSize(16, 16))
        self.addToolBar(toolbar)

        submenu_action = QAction(QIcon("document-import.png"), "Import", self)
        submenu_action.setStatusTip("Import file")
        submenu_action.triggered.connect(self.create_import_submenu)
        toolbar.addAction(submenu_action)
        toolbar.addSeparator()

        toolbar.addAction(self.create_export_action())
        toolbar.addSeparator()

        toolbar.addAction(self.show_history_action())

        self.setStatusBar(QStatusBar(self))

    def setup_menu(self):
        menu = self.menuBar()
        file_menu = menu.addMenu("&File")

        file_submenu = file_menu.addMenu("Import")
        file_submenu.addAction(self.create_import_action_megan())
        file_submenu.addAction(self.create_import_action_other())
        file_submenu.addAction(self.create_import_action_metadata())
        file_menu.addSeparator()

        file_menu.addAction(self.create_export_action())
        file_menu.addSeparator()

        file_menu.addAction(self.import_from_server_action())
        file_menu.addSeparator()

        file_submenu = file_menu.addMenu("Submenu")
        file_submenu.addAction(self.show_history_action())

        file_menu = menu.addMenu("&Window")

        file_menu = menu.addMenu("&Help")
        file_menu.addAction(self.create_help_action())
        file_menu.addSeparator()

    def create_import_submenu(self):
        menu = QMenu(self)

        action1 = QAction("MEGAN table", self)
        action1.triggered.connect(lambda: self.open_file_dialog())
        menu.addAction(action1)

        action2 = QAction("Other table", self)
        action2.triggered.connect(lambda: self.open_file_dialog())
        menu.addAction(action2)

        metadata = QAction("Metadata", self)
        metadata.triggered.connect(lambda: self.open_file_dialog())
        menu.addAction(metadata)

        cursor_pos = self.mapFromGlobal(self.cursor().pos())
        menu.exec(self.mapToGlobal(cursor_pos))

    # Actions on the toolbar and menu
    def create_import_action_megan(self):
        import_action = QAction(QIcon("globe.png"), "MEGAN table", self)
        import_action.setStatusTip("MEGAN table")
        # You can enter keyboard shortcuts using key names (e.g. Ctrl+p)
        # Qt.namespace identifiers (e.g. Qt.CTRL + Qt.Key_P)
        # or system agnostic identifiers (e.g. QKeySequence.StandardKey.Print)

        import_action.triggered.connect(self.open_file_dialog)
        return import_action

    def create_import_action_other(self):
        import_action = QAction(QIcon("document-import.png"), "Other", self)
        import_action.setStatusTip("other")
        # You can enter keyboard shortcuts using key names (e.g. Ctrl+p)
        # Qt.namespace identifiers (e.g. Qt.CTRL + Qt.Key_P)
        # or system agnostic identifiers (e.g. QKeySequence.StandardKey.Print)

        import_action.triggered.connect(self.open_file_dialog)
        return import_action

    def create_import_action_metadata(self):
        import_action = QAction(QIcon("document-import.png"), "Import metadata", self)
        import_action.setStatusTip("Import metadata")
        # You can enter keyboard shortcuts using key names (e.g. Ctrl+p)
        # Qt.namespace identifiers (e.g. Qt.CTRL + Qt.Key_P)
        # or system agnostic identifiers (e.g. QKeySequence.StandardKey.Print)

        import_action.triggered.connect(self.open_file_dialog)
        return import_action

    def create_export_action(self):
        export_action = QAction(QIcon("inbox-upload.png"), "Export file", self)
        export_action.setStatusTip("Export file")
        export_action.triggered.connect(self.onMyToolBarButtonClick)
        return export_action

    def show_history_action(self):
        show_history_action = QAction(QIcon("clock.png"), "Show history", self)
        show_history_action.setStatusTip("Show history of uploads")
        show_history_action.triggered.connect(self.onMyToolBarButtonClick)
        return show_history_action

    def import_from_server_action(self):
        import_from_server_action = QAction(QIcon("server--arrow.png"), "Import from server", self)
        import_from_server_action.setStatusTip("Select and import from server")
        return import_from_server_action

    def create_help_action(self):
        help_action = QAction(QIcon("globe.png"), "Ask in forum", self)
        help_action.setStatusTip("Ask for help in forum!")
        return help_action

    # Functions triggered by actions
    def onMyToolBarButtonClick(self, s):
        print("click", s)

    def button_clicked(self, s):
        print("click", s)
        dlg = CustomDialog()
        dlg.exec()

    def open_file_dialog(self):
        dialog = QFileDialog(self)
        dialog.setDirectory(r'C:\Users\egoro\Downloads')
        dialog.setFileMode(QFileDialog.FileMode.ExistingFile)
        dialog.setNameFilter("Tab-sep. files (*.csv *.tsv)")
        dialog.setViewMode(QFileDialog.ViewMode.List)
        if dialog.exec():
            file_name = dialog.selectedFiles()
            if file_name:
                self.file_path = file_name[0]
                otu_mat, tax_mat, metadata = self.process_file(self.file_path)
                input = MicrobiomeDataAnalyzer(otu_mat, tax_mat, metadata)
                self.plot_stat(input)

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

    def plot_stat(self, input):

        DT = dt.Frame(input.OTU_table)
        self.label.setText(f"Selected file: {self.file_path}\nOTU Table: {DT}")

        # Create the table view
        self.table_view = QTableView()

        # Create a model
        model = QStandardItemModel(DT.nrows, DT.ncols)

        # Set headers
        model.setHorizontalHeaderLabels(DT.names)

        # Fill the model with data
        for i in range(DT.nrows):
            for j in range(DT.ncols):
                item = QStandardItem(str(DT[i, j]))
                model.setItem(i, j, item)

        # Assign the model to the table view
        self.table_view.setModel(model)

        # Optional: Set alternating row colors for better readability
        self.table_view.setAlternatingRowColors(True)

        sc_top = MplCanvas(self, dpi=100)
        # sc_rank = MplCanvas(self, width=3, height=2, dpi=100)
        # sc_pcoa = MplCanvas(self, width=3, height=2, dpi=100)

        # Set maximum size
        sc_top.setMaximumSize(800, 400)  # Set maximum width and height

        # Plot on canvas
        input.plot_top(5, sc_top)
        # input.plot_rank("Genus", sc_rank)
        # input.plot_pcoa(sc_pcoa)

        # Set layout
        layout = QVBoxLayout()
        layout.addWidget(sc_top)
        # layout.addWidget(sc_rank)
        # layout.addWidget(sc_pcoa)
        layout.addWidget(self.table_view)

        central_widget = QWidget()
        central_widget.setLayout(layout)

        self.setCentralWidget(central_widget)

    def open_single_file_dialog(self):
        filename, ok = QFileDialog.getOpenFileName(
            self,
            "Select a File",
            "D:\\icons\\avatar\\",
            "Images (*.png *.jpg)"
        )
        if filename:
            path = Path(filename)
            self.filename_edit.setText(str(path))

app = QApplication(sys.argv)
w = MainWindow()
w.show()
app.exec()

# Todo
# Import file and apply functions from thesis
# Separete view and model
# Separate functions in separate files and classes
# Create a backup (github)

# SUbmenu for import: megan or other
# One class for functions
# Include metadata

# Data tabels in python
# Figure canvas

# metadata: sample id, groups
# Qiime to metadata

#datatable
