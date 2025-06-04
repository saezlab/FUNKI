import tkinter as tk
from tkinter import ttk

from funki import __version__

from tabs import all_tabs
from style import load_style


class Funki:

    def __init__(self, root):

        root.title('FUNKI v%s' % __version__)

        # Main frame
        mainframe = ttk.Frame(root)
        mainframe.pack(fill='both', expand=True)

        # Tab manager
        tab_manager = ttk.Notebook(
            mainframe,
            width=500,
            height=800,
        )
        tab_manager.pack(fill='both', expand=True, pady=(50, 10))

        # Adding tabs
        tabs = {}

        for name, (n, tab) in all_tabs.items():

            tabs[n] = tab(tab_manager)
            
            for child in tabs[n].winfo_children(): 
            
                child.grid_configure(padx=5, pady=5)
            
            tab_manager.add(tabs[n], text=name)


if __name__ == '__main__':

    root = tk.Tk()
    style = load_style()

    app = Funki(root)
    root.mainloop()