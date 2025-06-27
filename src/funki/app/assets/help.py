import tkinter as tk
from tkinter import ttk

from funki import __version__

from funki.app.utils import ScrollText
from funki.app.utils import PATH_LOGO
from funki.app.assets.msg_help import HELPMSG


class Help(tk.Toplevel):

    def __init__(self):

        super().__init__()

        self.title('FUNKI manual')
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        self.mainframe = ttk.Frame(self)
        self.mainframe.grid(
            column=0,
            row=0,
            sticky='NSEW',
            padx=(15, 15),
            pady=(15, 15)
        )
        self.mainframe.columnconfigure(0, weight=1)
        self.mainframe.columnconfigure(1, weight=3)
        self.mainframe.rowconfigure(0, weight=0)
        self.mainframe.rowconfigure(1, weight=1)

        # TOC
        self.tocvals = tk.StringVar(value=list(HELPMSG.keys()))
        self.toc = tk.Listbox(
            self.mainframe,
            listvariable=self.tocvals,
            selectmode='single'
        )
        self.toc.grid(row=0, column=0, sticky='NSEW', rowspan=2)
        self.toc.bind('<<ListboxSelect>>', self.update)

        # Title
        self.title = ttk.Label(
            self.mainframe,
            anchor='nw',
            padding=(10, 0, 0, 0),
            font=('Arial', 18, 'bold'),
        )
        self.title.grid(row=0, column=1, sticky='NSEW')

        # Contents
        self.content = ScrollText(
            self.mainframe,
            width=50,
        )
        self.content.grid(row=1, column=1, sticky='NSEW')

        # Initial position
        self.toc.selection_set(0)
        self.update()


    def update(self, *ev):

        k = self.toc.get(self.toc.curselection()[0])
        self.title['text'] = k
        self.content.update(HELPMSG[k])


class About(tk.Toplevel):

    def __init__(self):

        super().__init__()
        self.title('About FUNKI')
        self.resizable(False, False)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        mainframe = ttk.Frame(self)
        mainframe.grid(
            column=0,
            row=0,
            sticky='NSEW',
            padx=(15, 15),
            pady=(15, 15)
        )
        mainframe.columnconfigure(0, weight=1)
        mainframe.rowconfigure(0, weight=0) # Header (logo)
        mainframe.rowconfigure(1, weight=0) # Title (name & version)
        mainframe.rowconfigure(2, weight=1) # Desctiption
        mainframe.rowconfigure(3, weight=0) # Close button

        logo = tk.PhotoImage(file=PATH_LOGO)
        header = ttk.Label(
            mainframe,
            image=logo,
            anchor='center',
            padding=(10, 10, 10, 10),
        )
        header.image = logo # Avoiding garbage collection
        header.grid(column=0, row=0, sticky='NSEW')

        title = ttk.Label(
            mainframe,
            anchor='center',
            font=('Arial', 18, 'bold'),
            text='FUNKI v%s' % __version__,
        )
        title.grid(row=1, column=0, sticky='NSEW')

        descript = ttk.Label(
            mainframe,
            anchor='center',
            font='Arial',
            justify='center',
            text=(
                'Developed by Nicolàs Palacio-Escat\n'
                'Based on previous work from\n'
                'Rosa Hernansaiz-Ballesteros and Hanna Schumacher\n\n'
                'Developed at Saezlab\n'
                'Institute for Computational Biomedicine\n'
                'University Hospital Heidelberg\n\n'
                'https://github.com/saezlab/FUNKI\n\n'
                'Licensed under GPL-3.0'
            ),
            padding=(0, 10, 0, 10)
        )
        descript.grid(row=2, column=0, sticky='NSEW')

        button_close = ttk.Button(
            mainframe,
            text='Close',
            command=self.destroy,
            padding=(5, 5, 5, 5),
        )
        button_close.grid(row=3, column=0)