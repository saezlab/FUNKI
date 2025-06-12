import tkinter as tk
from tkinter import ttk

from utils import read_text
from utils import WrapLabel


PATH_MSG = 'src/funki/app/assets/msg_home.txt'


class TabHome(ttk.Frame):

    def __init__ (self, parent, controller, **options):

        super().__init__(parent, **options)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=0)

        # Welcome text
        maintext = WrapLabel(
            self,
            text=read_text(PATH_MSG),
            anchor='nw',
        )
        maintext.grid(row=0, column=0, sticky='NSEW', columnspan=2)

        # Button frames
        lframe = ttk.Frame(self)
        lframe.grid(row=1, column=0)
        lframe.rowconfigure(0, weight=0)
        lframe.rowconfigure(1, weight=0)
        lframe.columnconfigure(0, weight=0)

        rframe = ttk.Frame(self)
        rframe.grid(row=1, column=1)
        rframe.rowconfigure(0, weight=0)
        rframe.rowconfigure(1, weight=0)
        rframe.columnconfigure(0, weight=0)

        # Buttons raw data
        button_loadraw = ttk.Button(
            lframe,
            text='Load data',
            padding=(5, 5, 5, 5),
        )
        button_loadraw.grid(row=0, column=0, pady=(5, 5))
        button_viewraw = ttk.Button(
            lframe,
            text='View data',
            padding=(5, 5, 5, 5),
        )
        button_viewraw.grid(row=1, column=0, pady=(5, 5))

        # Buttons metadata
        button_loadobs = ttk.Button(
            rframe,
            text='Load metadata',
            padding=(5, 5, 5, 5),
        )
        button_loadobs.grid(row=0, column=1, pady=(5, 5))
        button_viewobs = ttk.Button(
            rframe,
            text='View metadata',
            padding=(5, 5, 5, 5),
        )
        button_viewobs.grid(row=1, column=1, pady=(5, 5))