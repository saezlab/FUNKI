import os

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog as fd
import tksvg
import pandas as pd

from funki import __version__
from funki.input import read

from utils import PATH_LOGO
from tabs import TABS
from style import load_style
from assets.help import Help, About


class FunkiApp:

    _platform = None

    def __init__(self, root):
        '''
        Establishes application's root and platform, runs all setup methods and
        starts a new project.
        '''

        self._root = root
        self._platform = self._root.tk.call('tk', 'windowingsystem')

        self.new_project()


    def _setup_root(self):
        '''
        Sets up the application's root configuration.
        '''

        self._root.title('FUNKI v%s' % __version__)
        self._root.option_add('*tearOff', False) # Avoids menu detachment
        self._root.columnconfigure(0, weight=1)
        self._root.rowconfigure(0, weight=1)


    def _setup_menu(self):
        '''
        Sets up the top level menus.
        '''

        self.menubar = tk.Menu(self._root)
        self._root.config(menu=self.menubar)

        self.menu_file = tk.Menu(self.menubar)
        self.menubar.add_cascade(menu=self.menu_file, label='File')
        self.menu_help = tk.Menu(self.menubar)
        self.menubar.add_cascade(menu=self.menu_help, label='Help')

        menu_options = [
            (
                self.menu_file,
                'New project',
                self.new_project
            ),
            (
                self.menu_file,
                'Load data',
                lambda: self.open_file(dtype='raw')
            ),
            (
                self.menu_file,
                'Load metadata',
                lambda: self.open_file(dtype='obs')
            ),
            (
                self.menu_help,
                'FUNKI manual',
                self.open_manual
            ),
            (
                self.menu_help,
                'About FUNKI',
                self.open_about
            ),
        ]

        for menu, label, command in menu_options:

                menu.add_command(label=label, command=command)

        # Disable metadata loading until data is available
        self.menu_file.entryconfig('Load metadata', state='disabled')


    def _setup_mainframe(self):
        '''
        Sets up the application's main parent frame where all children widgets
        are to be placed. Calls setup of main children: header and notebook.
        '''

        self.mainframe = ttk.Frame(self._root)
        self.mainframe.grid(
            column=0,
            row=0,
            sticky='NSEW',
            padx=(5, 5),
            pady=(5, 5)
        )
        self.mainframe.columnconfigure(0, weight=1)
        self.mainframe.rowconfigure(0, weight=0)
        self.mainframe.rowconfigure(1, weight=1)

        self._setup_header()
        self._setup_notebook()


    def _setup_header(self):
        '''
        Sets up the header logo.
        '''

        logo = tksvg.SvgImage(file=PATH_LOGO, scale=0.25)

        header = ttk.Label(self.mainframe, image=logo, padding=(10, 10, 10, 10))
        header.image = logo # Avoiding garbage collection
        header.grid(column=0, row=0, sticky='NSEW')


    def _setup_notebook(self):
        '''
        Sets up the notebook, which is the main window of the application, where
        all the tabs for the different steps are shown. These are coded in the
        different submodules on the `tabs/` folder.
        '''

        self.tab_manager = ttk.Notebook(
            self.mainframe,
            width=500,
            height=800,
        )
        self.tab_manager.grid(column=0, row=1, sticky='NSEW')

        # Adding tabs
        tabs = {}

        for name, (n, tab) in TABS.items():

            tabs[n] = tab(self.tab_manager)
            
            for child in tabs[n].winfo_children():
            
                child.grid_configure(padx=5, pady=5)
            
            self.tab_manager.add(tabs[n], text=name, sticky='NSEW')


    def new_project(self):
        '''
        Clears current data and resets the application.
        '''

        self.data = None

        self._setup_root()
        self._setup_menu()
        self._setup_mainframe()


    def open_file(self, dtype=None):
        '''
        Pops up a open file dialog and passes the path to FUNKI's read method.
        '''

        path = fd.askopenfilename()

        # Loading of measurement data
        if path and dtype == 'raw':

            self.data = read(path)

        # Loading of metadata
        elif path and dtype == 'obs' and self.data:

            fname, ext = os.path.splitext(path)

            self.data.obs = pd.read_csv(
                path,
                sep=',' if ext == '.csv' else '\t',
                index_col=0
            )

        # Reactivate the Load metadata in menu
        if self.data and dtype == 'raw':

            self.menu_file.entryconfig('Load metadata', state='normal')


    def open_manual(self):

        Help()


    def open_about(self):

        pass


if __name__ == '__main__':

    root = tk.Tk()
    load_style(root)

    app = FunkiApp(root)
    root.mainloop()