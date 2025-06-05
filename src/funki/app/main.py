import tkinter as tk
from tkinter import ttk
import tksvg

from funki import __version__

from tabs import TABS
from style import load_style


PATH_LOGO = '../assets/logos/funki_logo.svg'

class Funki:

    _root = None
    _platform = None

    def __init__(self, root):
        self._root = root

        self._setup_root()
        self._setup_menu()
        self._setup_mainframe()


    def _setup_root(self):

        self._platform = self._root.tk.call('tk', 'windowingsystem')

        self._root.title('FUNKI v%s' % __version__)
        self._root.option_add('*tearOff', False) # Avoids menu detachment
        self._root.columnconfigure(0, weight=1)
        self._root.rowconfigure(0, weight=1)


    def _setup_menu(self):

        menubar = tk.Menu(self._root)
        self._root.config(menu=menubar)

        menu_file = tk.Menu(menubar)
        menubar.add_cascade(menu=menu_file, label='File')

        menu_file.add_command(label='Open...', command=self.openFile)


    def openFile(self):
        pass


    def _setup_mainframe(self):

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

        logo = tksvg.SvgImage(file=PATH_LOGO, scale=0.25)

        header = ttk.Label(self.mainframe, image=logo, padding=(10, 10, 10, 10))
        header.image = logo # Avoiding garbage collection
        header.grid(column=0, row=0, sticky='NSEW')


    def _setup_notebook(self):

        tab_manager = ttk.Notebook(
            self.mainframe,
            width=500,
            height=800,
        )
        tab_manager.grid(column=0, row=1, sticky='NSEW')

        # Adding tabs
        tabs = {}

        for name, (n, tab) in TABS.items():

            tabs[n] = tab(tab_manager)
            
            for child in tabs[n].winfo_children():
            
                child.grid_configure(padx=5, pady=5)
            
            tab_manager.add(tabs[n], text=name, sticky='NSEW')


if __name__ == '__main__':

    root = tk.Tk()
    load_style(root)

    app = Funki(root)
    root.mainloop()