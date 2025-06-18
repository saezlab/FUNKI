import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt

from utils import Figure
from utils import  LabeledWidget


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
        title_dex.grid(row=0, columnspan=2, sticky='NSWE')

        # Contrast variable selector
        self.combox_contrast_var = LabeledWidget(
            self,
            ttk.Combobox,
            'Select group variable for the contrast: ',
            lpos='n',
            wget_kwargs={'state': 'disabled'},
            wget_grid_kwargs={'sticky': 'EW', 'weight': 1},
            label_grid_kwargs={'sticky': 'EW', 'weight': 0},
        )
        self.combox_contrast_var.grid(row=1, columnspan=2, sticky='NSEW')
        self.combox_contrast_var.wg.bind('<<ComboboxSelected>>', self._update)

        # Contrast group A
        self.combox_A = LabeledWidget(
            self,
            ttk.Combobox,
            'Group of samples to contrast from: ',
            lpos='n',
            wget_kwargs={'state': 'disabled'},
            wget_grid_kwargs={'sticky': 'EW', 'weight': 1},
            label_grid_kwargs={'sticky': 'EW', 'weight': 0},
        )
        self.combox_A.grid(row=2, column=0, sticky='NSWE')
        self.combox_A.wg.bind('<<ComboboxSelected>>', self._update)

        # Contrast group B
        self.combox_B = LabeledWidget(
            self,
            ttk.Combobox,
            'Group of samples to contrast against: ',
            lpos='n',
            wget_kwargs={'state': 'disabled'},
            wget_grid_kwargs={'sticky': 'EW', 'weight': 1},
            label_grid_kwargs={'sticky': 'EW', 'weight': 0},
        )
        self.combox_B.grid(row=2, column=1, sticky='NSWE')
        self.combox_B.wg.bind('<<ComboboxSelected>>', self._update)


    def _update(self, *ev):

        if self.controller.data and not self.controller.data.obs.empty:

            obs_key = self.combox_contrast_var.wg.get()

            # Main combobox is empty
            if not obs_key:

                obs_keys = sorted([
                    c for c in self.controller.data.obs_keys()
                    if all([
                        isinstance(i, str)
                        for i in self.controller.data.obs[c]
                    ])
                ])

                if obs_keys:

                    obs_key = obs_keys[0]
                    self.combox_contrast_var.wg.configure(
                        state='readonly',
                        values=obs_keys,
                    )
                    self.combox_contrast_var.wg.set(obs_key)

            # Main combobox not empty
            if obs_key:

                options = sorted(set(self.controller.data.obs[obs_key]))

                # A/B are empty or new obs_key
                if (
                    (not self.combox_A.wg.get() or not self.combox_B.wg.get())
                    or (
                        self.combox_A.wg.get() not in options
                        or self.combox_A.wg.get() not in options
                    )
                ):

                    self.combox_A.wg.configure(
                        state='readonly',
                        values=[n for i, n in enumerate(options) if i != 1],
                    )
                    self.combox_A.wg.set(options[0])

                    self.combox_B.wg.configure(
                        state='readonly',
                        values=[n for i, n in enumerate(options) if i != 0],
                    )
                    self.combox_B.wg.set(options[1])

                # A/B not empty -> ensure mutual exclusivity
                else:
                    curA = self.combox_A.wg.get()
                    curB = self.combox_B.wg.get()

                    self.combox_A.wg.configure(
                        values=[i for i in options if i != curB],
                    )

                    self.combox_B.wg.configure(
                        values=[i for i in options if i != curA],
                    )
