import tkinter as tk
from tkinter import ttk

from funki.app.utils import read_text
from funki.app.utils import WrapLabel
from funki.app.utils import PATH_MSG


class TabHome(ttk.Frame):

    def __init__ (self, parent, controller, **options):

        super().__init__(parent, **options)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=0)

        self.controller = controller

        # Welcome text
        WrapLabel(
            self,
            text=read_text(PATH_MSG),
            anchor='nw',
        ).grid(row=0, column=0, sticky='NSEW', columnspan=2)

        # Button frames
        lframe = ttk.Frame(self)
        lframe.grid(row=1, column=0)
        lframe.rowconfigure(0, weight=1)
        lframe.rowconfigure(1, weight=0)
        lframe.rowconfigure(2, weight=0)
        lframe.columnconfigure(0, weight=1)

        rframe = ttk.Frame(self)
        rframe.grid(row=1, column=1)
        rframe.rowconfigure(0, weight=1)
        rframe.rowconfigure(1, weight=0)
        rframe.rowconfigure(2, weight=0)
        rframe.columnconfigure(0, weight=1)

        # Buttons raw data
        ttk.Label(
            lframe,
            text='Raw data',
            style='Title.TLabel',
        ).grid(row=0, column=0, sticky='NSEW')
        self.button_loadraw = ttk.Button(
            lframe,
            text='Load',
            padding=(5, 5, 5, 5),
            command=lambda: self.controller.open_file(dtype='raw'),
        )
        self.button_loadraw.grid(row=1, column=0, pady=(5, 5))
        self.button_viewraw = ttk.Button(
            lframe,
            text='View',
            padding=(5, 5, 5, 5),
            command=lambda: self.controller.view_data(dtype='raw'),
            state='disabled',
        )
        self.button_viewraw.grid(row=2, column=0, pady=(5, 5))

        # Buttons metadata
        ttk.Label(
            rframe,
            text='Metadata',
            style='Title.TLabel',
        ).grid(row=0, column=0, sticky='NW')
        self.button_loadobs = ttk.Button(
            rframe,
            text='Load',
            padding=(5, 5, 5, 5),
            command=lambda: self.controller.open_file(dtype='obs'),
            state='disabled',
        )
        self.button_loadobs.grid(row=1, column=0, pady=(5, 5))
        self.button_viewobs = ttk.Button(
            rframe,
            text='View',
            padding=(5, 5, 5, 5),
            command=lambda: self.controller.view_data(dtype='obs'),
            state='disabled',
        )
        self.button_viewobs.grid(row=2, column=0, pady=(5, 5))


    def _update(self):

        if self.controller.data:

            self.button_viewraw.configure(state='normal')
            self.button_loadobs.configure(state='normal')

            if not self.controller.data.obs.empty:

                self.button_viewobs.configure(state='normal')
