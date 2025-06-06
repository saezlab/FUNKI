import tkinter as tk
from tkinter import ttk

from utils import WrapLabel
from .msg_help import HELPMSG


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
        self.content = WrapLabel(
            self.mainframe,
            width=50,
            anchor='nw',
            padding=(10, 0, 0, 0),
        )
        self.content.grid(row=1, column=1, sticky='NSEW')

        # Initial position
        self.toc.selection_set(0)
        self.update()


    def update(self, *ev):

        k = self.toc.get(self.toc.curselection()[0])
        self.title['text'] = k
        self.content['text'] = HELPMSG[k]