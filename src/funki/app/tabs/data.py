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
            table_data = Table(self, df, padding=(5, 5, 5, 5))
            table_data.grid(row=0, column=0, sticky='NSWE')

            if not self.controller.data.obs.empty:

                obs = self.controller.data.obs
                table_obs = Table(self, obs, padding=(5, 5, 5, 5))
                table_obs.grid(row=0, column=1, sticky='NSWE')
