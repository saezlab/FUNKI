import tkinter as tk
from tkinter import ttk

from utils import Table

import pandas as pd
import numpy as np


class TabData(ttk.Frame):

    def __init__ (self, parent, **options):

        super().__init__(parent, **options)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)

        df = pd.DataFrame(np.arange(12).reshape(3, 4), columns=list('ABCD'))

        table_data = Table(self, df)
        table_data.grid(row=0, column=0, sticky='NSEW')

        