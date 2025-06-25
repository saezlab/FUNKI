import os

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog as fd
import tksvg
import pandas as pd

from funki import __version__
from funki.input import read

from utils import PATH_LOGO
from utils import PopUpTable
from utils import check_num
from tabs import TABS
from style import load_style
from assets.help import Help
from assets.help import About


class FunkiApp(tk.Tk):
    # TODO: Check if some of the Label can be replaced with Text for better
    # formatting and handling

    _platform = None

    def __init__(self):
        '''
        Establishes application's root and platform, runs all setup methods and
        starts a new project.
        '''

        super().__init__()

        self.geometry('800x800')
        self._platform = self.tk.call('tk', 'windowingsystem')

        self.check_num = (self.register(check_num), '%P')

        self.new_project()


    def _setup_root(self):
        '''
        Sets up the application's root configuration.
        '''

        self.title('FUNKI v%s' % __version__)
        self.option_add('*tearOff', False) # Avoids menu detachment
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)


    def _setup_menu(self):
        '''
        Sets up the top level menus.
        '''

        self.menubar = tk.Menu(self)
        self.config(menu=self.menubar)

        self.menu_file = tk.Menu(self.menubar)
        self.menu_open = tk.Menu(self.menu_file)
        self.menu_save = tk.Menu(self.menu_file)
        self.menu_view = tk.Menu(self.menubar)
        self.menu_help = tk.Menu(self.menubar)

        # Keys are parent menu instance, values are either:
        # * '---': Will add a separator in that position
        # * tuple: Contains 3 elements:
        #   - Label/text, status and command will create that option
        #   - Label/text, None, parent menu instance will add a cascade
        menu_options = {
            self.menubar: [
                ('File', None, self.menu_file),
                ('View', None, self.menu_view),
                ('Help', None, self.menu_help),
            ],
            self.menu_file: [
                (
                    'New project',
                    'normal',
                    self.new_project
                ),
                ('Open...', None, self.menu_open),
                ('Save...', None, self.menu_save),
                '---',
                (
                    'Exit',
                    'normal',
                    self.destroy
                ),
            ],
            self.menu_open: [
                (
                    'Data',
                    'normal',
                    lambda: self.open_file(dtype='raw')
                ),
                (
                    'Metadata',
                    'disabled',
                    lambda: self.open_file(dtype='obs')
                ),
            ],
            self.menu_save: [
                (
                    'Data',
                    'disabled',
                    lambda: self.save_file(dtype='raw')
                ),
                (
                    'Metadata',
                    'disabled',
                    lambda: self.save_file(dtype='obs')
                ),
                (
                    'Differential expression',
                    'disabled',
                    lambda: self.save_file(dtype='dex')
                ),
                '---',
                (
                    'Configuration',
                    'disabled',
                    lambda: self.save_file(dtype='config')
                ),
            ],
            self.menu_view: [
                (
                    'Data',
                    'disabled',
                    lambda: self.view_data(dtype='raw')
                ),
                (
                    'Metadata',
                    'disabled',
                    lambda: self.view_data(dtype='obs')
                ),
                (
                    'Differential expression',
                    'disabled',
                    lambda: self.view_data(dtype='dex')
                ),
                (
                    'Gene Set collection',
                    'normal',
                    lambda: self.view_data(dtype='gsc')
                ),
            ],
            self.menu_help: [
                (
                    'FUNKI manual',
                    'normal',
                    self.open_manual
                ),
                (
                    'About FUNKI',
                    'normal',
                    self.open_about
                ),
            ],
        }

        for menu, options in menu_options.items():

            for label, state, command in options:

                if label == '-':

                    menu.add_separator()

                elif state == None:

                    menu.add_cascade(menu=command, label=label)

                else:

                    menu.add_command(label=label, command=command, state=state)


    def _setup_mainframe(self):
        '''
        Sets up the application's main parent frame where all children widgets
        are to be placed. Calls setup of main children: header and notebook.
        '''

        self.mainframe = ttk.Frame(self)
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

        self.tab_manager = ttk.Notebook(self.mainframe)
        self.tab_manager.grid(column=0, row=1, sticky='NSEW')

        # Adding tabs
        self.tabs = {}

        for name, (n, tab) in TABS.items():

            self.tabs[n] = tab(self.tab_manager, self)

            for child in self.tabs[n].winfo_children():

                child.grid_configure(padx=10, pady=10)

            self.tab_manager.add(self.tabs[n], text=name, sticky='NSEW')


    def _update(self):
        '''
        Handles the activation of diverse options upon loading of raw/obs data.
        '''

        if self.data:

            self.menu_open.entryconfig('Metadata', state='normal')
            self.menu_view.entryconfig('Data', state='normal')
            self.menu_save.entryconfig('Data', state='normal')
            self.menu_save.entryconfig('Configuration', state='normal')

            if not self.data.obs.empty:

                self.menu_view.entryconfig('Metadata', state='normal')
                self.menu_save.entryconfig('Metadata', state='normal')

        for tab in self.tabs.values():

            if hasattr(tab, '_update'):

                tab._update()


    def new_project(self):
        '''
        Clears current data and resets the application.
        '''

        self.data = None
        self.fname_raw = None
        self.fname_obs = None

        self._setup_root()
        self._setup_menu()
        self._setup_mainframe()


    def open_file(self, dtype=None):
        '''
        Pops up a open file dialog and passes the path to FUNKI's read method.
        '''

        filetypes = [
            ('All files', '*.*'),
            ('CSV files', '*.csv'),
            ('TSV files', '*.tsv'),
            ('TXT files', '*.txt'),
        ]
        filetypes += [
            ('H5AD files', '*.h5ad'),
            ('H5 files', '*.h5'),
            ('Excel files', '*.xlsx'),
            ('Loom files', '*.loom'),
            ('MTX files', '*.mtx'),
            ('UMI files', '*.gz'),
        ] if dtype == 'raw' else []

        path = fd.askopenfilename(defaultextension='.*', filetypes=filetypes)

        if not path:

            return

        # Loading of measurement data
        elif dtype == 'raw': # TODO: Make a copy in raw layer?

            self.fname_raw = path
            self.data = read(path)

        # Loading of metadata
        elif dtype == 'obs' and self.data:

            fname, ext = os.path.splitext(path)

            self.fname_obs = path
            obs = pd.read_csv(
                path,
                sep=',' if ext == '.csv' else '\t',
                index_col=0
            )
            # TODO: Handle samples not in data
            self.data.obs = self.data.obs.merge(
                obs,
                how='outer',
                left_index=True,
                right_index=True,
            )

        self._update()


    def save_file(self, dtype=None):
        '''
        Pops up a save file dialog and saves the corresponding file in the path.
        '''

        defaultextension = '.json' if dtype == 'config' else '.csv'
        filetypes = [
            ('JSON files', '*.json'),
        ] if dtype == 'config' else [
            ('CSV files', '*.csv'),
            ('TSV files', '*.tsv'),
            ('TXT files', '*.txt'),
        ]

        path = fd.asksaveasfilename(
            defaultextension=defaultextension,
            filetypes=filetypes + [('All files', '*.*')]
        )

        if not path:

            return

        fname, ext = os.path.splitext(path)
        sep = ',' if ext == '.csv' else '\t'

        if dtype == 'raw':

            df = self.data.to_df()

        elif dtype == 'obs':

            df = self.data.obs

        elif dtype == 'dex':

            df = self.data.var[[
                'baseMean',
                'log2FoldChange',
                'stat',
                'pvalue',
                'padj'
            ]]

        # Storing file
        if dtype == 'config':

            self.data.save_params(path)

        else:

            df.to_csv(path, sep=sep)


    def open_manual(self):

        Help()


    def open_about(self):

        About()


    def view_data(self, dtype=None):

        df = pd.DataFrame()
        title = ''

        if dtype == 'raw' and self.data:

            title = f'Data - {self.fname_raw}'
            df = self.data.to_df()

        elif dtype == 'obs' and (self.data and not self.data.obs.empty):

            title = f'Metadata - {self.fname_obs}'
            df = self.data.obs

        elif dtype == 'dex' and 'diff_exp' in self.data.uns['funki']:

            a = self.tabs['dex'].groupA.get()
            b = self.tabs['dex'].groupB.get()

            title = f'DEX - {a} vs. {b}'
            df = self.data.varm[f'{a}_vs_{b}'] # TODO: add contingency measures

        elif dtype == 'gsc':

            if self.tabs['enrich'].net.empty:

                self.tabs['enrich']._get_resource()

            res = self.tabs['enrich'].gset.get()
            org = self.tabs['enrich'].org.get()

            df = self.tabs['enrich'].net
            title = f'{res} ({org})'

        command = (
            (lambda: self.save_file(dtype=dtype))
            if dtype in {'raw', 'obs', 'dex'}
            else None
        )

        PopUpTable(df, title=title, save_command=command)


if __name__ == '__main__':

    app = FunkiApp()
    load_style(app)

    app.mainloop()