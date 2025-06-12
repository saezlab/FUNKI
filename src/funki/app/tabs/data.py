import tkinter as tk
from tkinter import ttk

import pandas as pd
import numpy as np


class TabData(ttk.Frame):

    def __init__ (self, parent, controller, **options):

        super().__init__(parent, **options)

        self.controller = controller

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
