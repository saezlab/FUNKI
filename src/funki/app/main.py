import os

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog as fd
import tksvg
import pandas as pd

from funki import __version__
from funki.input import read

from utils import PATH_LOGO, Table, check_num
from tabs import TABS
from style import load_style
from assets.help import Help, About


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

        self.geometry('700x800')
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
                '---',
                (
                    'Exit',
                    'normal',
                    self.destroy
                ),
            ],
            self.menu_open: [
                (
                    'Load data',
                    'normal',
                    lambda: self.open_file(dtype='raw')
                ),
                (
                    'Load metadata',
                    'disabled',
                    lambda: self.open_file(dtype='obs')
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


    def _update(self, dtype=None):
        '''
        Handles the activation of diverse options upon loading of raw/obs data.
        '''

        if dtype == 'raw' and self.data:

            self.menu_open.entryconfig('Load metadata', state='normal')
            self.menu_view.entryconfig('Data', state='normal')

        elif dtype == 'obs' and (self.data and not self.data.obs.empty):

            self.menu_view.entryconfig('Metadata', state='normal')

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

        path = fd.askopenfilename()

        # Loading of measurement data
        if path and dtype == 'raw':

            self.fname_raw = path
            self.data = read(path)

        # Loading of metadata
        elif path and dtype == 'obs' and self.data:

            fname, ext = os.path.splitext(path)

            self.fname_obs = path
            self.data.obs = pd.read_csv(
                path,
                sep=',' if ext == '.csv' else '\t',
                index_col=0
            )

        self._update(dtype=dtype)


    def open_manual(self):

        Help()


    def open_about(self):

        About()


    def view_data(self, dtype=None): # NOTE: Could move to new class?

        df = pd.DataFrame()
        title = ''

        if dtype == 'raw' and self.data:

            title = 'Data - %s' % self.fname_raw
            df = self.data.to_df()

        elif dtype == 'obs' and (self.data and not self.data.obs.empty):

            title = 'Metadata - %s' % self.fname_obs
            df = self.data.obs

        win = tk.Toplevel()

        win.title(title)
        win.columnconfigure(0, weight=1)
        win.rowconfigure(0, weight=1)
        win.rowconfigure(1, weight=0)

        if not df.empty:

            table = Table(win, df, padding=(5, 5, 5, 5))
            table.grid(row=0, column=0, sticky='NSWE')

        button_close = ttk.Button(
            win,
            text='Close',
            command=win.destroy,
            padding=(5, 5, 5, 5),
        )
        button_close.grid(row=1, column=0)


if __name__ == '__main__':

    app = FunkiApp()
    load_style(app)

    app.mainloop()