import tkinter as tk
from tkinter import ttk

from funki import __version__

from tabs import all_tabs

class Funki:

    def __init__(self, root):

        self.root = root
        self.root.title('FUNKI v%s' % __version__)

        self.tab_manager = ttk.Notebook(root, width=500, height=800)
        self.tab_manager.pack(expand=1, fill='both')
        self.tabs = {}

        for name, (n, tab) in all_tabs.items():

            self.tabs[n] = tab(self.tab_manager)
            self.tab_manager.add(self.tabs[n], text=name)


if __name__ == '__main__':

    root = tk.Tk()

    app = Funki(root)
    root.mainloop()