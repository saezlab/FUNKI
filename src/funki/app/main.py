import tkinter as tk
from tkinter import ttk

from funki import __version__

from tabs import all_tabs
from style import load_style


class Funki:

    def __init__(self, root):

        self.root = root
        self.root.title('FUNKI v%s' % __version__)

        # Main frame
        self.mainframe = ttk.Frame(root)
        self.mainframe.pack(fill='both', expand=True)

        # Tab manager
        self.tab_manager = ttk.Notebook(
            self.mainframe,
            width=500,
            height=800,
        )
        self.tab_manager.pack(fill='both', expand=True, pady=(50, 10))

        # Adding tabs
        self.tabs = {}

        for name, (n, tab) in all_tabs.items():

            self.tabs[n] = tab(self.tab_manager)
            
            for child in self.tabs[n].winfo_children(): 
            
                child.grid_configure(padx=5, pady=5)
            
            self.tab_manager.add(self.tabs[n], text=name)


if __name__ == '__main__':

    root = tk.Tk()
    style = load_style()

    app = Funki(root)
    root.mainloop()