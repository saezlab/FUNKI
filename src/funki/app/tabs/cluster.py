import tkinter as tk
from tkinter import ttk


class TabClust(ttk.Frame):

    def __init__ (self, parent, controller, **options):

        super().__init__(parent, **options)

        self.controller = controller