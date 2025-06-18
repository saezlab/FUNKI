import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt

from funki.plots import plot_obs

from utils import Figure
from utils import LabeledWidget


class TabData(ttk.Frame):

    def __init__ (self, parent, controller, **options):

        super().__init__(parent, **options)

        self.controller = controller

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=0)
        self.rowconfigure(1, weight=0)
        self.rowconfigure(2, weight=1)

        # Raw data panel
        title_raw = ttk.Label(
            self,
            text='Raw data:',
            style='Title.TLabel',
            width=20
        )
        title_raw.grid(row=0, column=0, sticky='NSWE')
        summary = ttk.Frame(self)
        summary.rowconfigure(0, weight=0)
        summary.rowconfigure(1, weight=0)

        self.nobs = tk.StringVar()
        self.nvar = tk.StringVar()

        LabeledWidget(
            summary,
            ttk.Label,
            'No. of observations: ',
            lpos='n',
            wget_kwargs={'textvariable': self.nobs},
            wget_grid_kwargs={'sticky': 'w', 'weight': 0},
            label_grid_kwargs={'sticky': 'w', 'weight': 0},
        ).grid(row=0, sticky='W')

        LabeledWidget(
            summary,
            ttk.Label,
            'No. of variables: ',
            lpos='n',
            wget_kwargs={'textvariable': self.nvar},
            wget_grid_kwargs={'sticky': 'w', 'weight': 0},
            label_grid_kwargs={'sticky': 'w', 'weight': 0},
        ).grid(row=1, sticky='W')

        summary.grid(column=0, row=2, sticky='NSWE')

        # Obs data panel
        title_obs = ttk.Label(
            self,
            text='Metadata:',
            style='Title.TLabel',
            width=20
        )
        title_obs.grid(row=0, column=1, sticky='NSWE')
        self.combox = LabeledWidget(
            self,
            ttk.Combobox,
            'Select variable to visualize: ',
            lpos='w',
            wget_kwargs={'state': 'disabled'},
            wget_grid_kwargs={'sticky': 'w', 'weight': 0},
            label_grid_kwargs={'sticky': 'w', 'weight': 0},
        )
        self.combox.grid(row=1, column=1, sticky='NSWE')
        self.combox.wg.bind('<<ComboboxSelected>>', self._update)

        self.fig, self.ax = plt.subplots()

        self.figframe = Figure(self, self.fig)
        self.figframe.grid(row=2, column=1, sticky='NSWE')


    def _update(self, *ev):

        if self.controller.data:

            no, nv = self.controller.data.shape
            self.nobs.set(str(no))
            self.nvar.set(str(nv))

            if not self.controller.data.obs.empty:

                # Set combobox
                if not self.combox.wg.get():

                    obs_keys = sorted(self.controller.data.obs_keys())

                    self.combox.wg.configure(
                        state='readonly',
                        values=obs_keys,
                    )
                    self.combox.wg.set(obs_keys[0])

                # Plot
                self.ax.clear()
                plot_obs(
                    self.controller.data,
                    obs_var=self.combox.wg.get(),
                    ax=self.ax
                )
                self.figframe._update()
