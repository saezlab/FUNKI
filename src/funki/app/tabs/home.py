import tkinter as tk
from tkinter import ttk

from utils import read_text
from utils import WrapLabel


PATH_MSG = 'src/funki/app/assets/msg_home.txt'


class TabHome(ttk.Frame):

    def __init__ (self, parent, **options):

        super().__init__(parent, **options)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        maintext = WrapLabel(
            self,
            text=read_text(PATH_MSG),
            anchor='nw',
        )
        maintext.grid(row=0, column=0, sticky='NSEW')
