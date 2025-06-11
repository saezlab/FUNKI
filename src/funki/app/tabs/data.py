import tkinter as tk
from tkinter import ttk

from utils import Table

import pandas as pd
import numpy as np


class TabData(ttk.Frame):

    def __init__ (self, parent, controller, **options):

        super().__init__(parent, **options)

        self.controller = controller

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)


    def update(self):

        if self.controller.data:

            df = self.controller.data.to_df()
            table_data = Table(self, df)
            table_data.grid(row=0, column=0)
