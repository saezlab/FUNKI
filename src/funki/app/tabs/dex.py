import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt

from utils import Figure, LabeledWidget


class TabDex(ttk.Frame):

    def __init__ (self, parent, controller, **options):

        super().__init__(parent, **options)

        self.controller = controller

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=0)
        self.rowconfigure(1, weight=0)
        self.rowconfigure(2, weight=0)
        self.rowconfigure(3, weight=1)

        title_dex = ttk.Label(
            self,
            text='Differential expression:',
            style='Title.TLabel',
            width=20
        )
        title_dex.grid(row=0, column=0, sticky='NSWE')

        self.combox_contrast_var = LabeledWidget(
            self,
            ttk.Combobox,
            'Select group variable for the contrast: ',
            lpos='w',
            wget_kwargs={'state': 'disabled'}
        )
        self.combox_contrast_var.grid(row=1, columnspan=2, sticky='NSWE')
        self.combox_contrast_var.wg.bind('<<ComboboxSelected>>', self._update)


    def _update(self, *ev):

        if self.controller.data and not self.controller.data.obs.empty:

            # Set combobox
            if not self.combox_contrast_var.wg.get():

                obs_keys = list(self.controller.data.obs_keys())

                self.combox_contrast_var.wg.configure(
                    state='readonly',
                    values=obs_keys,
                )
                self.combox_contrast_var.wg.set(obs_keys[0])