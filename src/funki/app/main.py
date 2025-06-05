import tkinter as tk
from tkinter import ttk
import tksvg

from funki import __version__

from tabs import all_tabs
from style import load_style


class Funki:

    def __init__(self, root):

        root.title('FUNKI v%s' % __version__)
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)

        # Main frame
        mainframe = ttk.Frame(root)
        mainframe.grid(column=0, row=0, sticky='NSEW')
        mainframe.columnconfigure(0, weight=1)
        mainframe.rowconfigure(0, weight=0)
        mainframe.rowconfigure(1, weight=1)

        logo = tksvg.SvgImage(
            file='../assets/logos/funki_logo.svg',
            scale=0.25
        )

        header = ttk.Label(mainframe, image=logo, padding=(10, 10, 10, 10))
        header.image = logo # Avoiding garbage collection
        header.grid(column=0, row=0, sticky='NSEW')

        # Tab manager
        tab_manager = ttk.Notebook(
            mainframe,
            width=500,
            height=800,
        )
        tab_manager.grid(column=0, row=1, sticky='NSEW')

        # Adding tabs
        tabs = {}

        for name, (n, tab) in all_tabs.items():

            tabs[n] = tab(tab_manager)
            
            for child in tabs[n].winfo_children():
            
                child.grid_configure(padx=5, pady=5)
            
            tab_manager.add(tabs[n], text=name, sticky='NSEW')


if __name__ == '__main__':

    root = tk.Tk()
    load_style(root)

    app = Funki(root)
    root.mainloop()