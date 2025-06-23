import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt

from funki import plots

from utils import Figure
from utils import LabeledWidget


class TabClust(ttk.Frame):

    def __init__ (self, parent, controller, **options):

        super().__init__(parent, **options)

        self.controller = controller

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=0)
        self.rowconfigure(1, weight=0)
        self.rowconfigure(2, weight=0)
        self.rowconfigure(3, weight=0)
        self.rowconfigure(4, weight=0)
        self.rowconfigure(5, weight=0)
        self.rowconfigure(6, weight=1)

        # Embedding panel
        ttk.Label(
            self,
            text='Embedding:',
            style='Title.TLabel',
            width=20
        ).grid(row=0, column=0, sticky='NSWE')

        self.harmony = tk.BooleanVar()
        LabeledWidget(
            self,
            ttk.Checkbutton,
            'Apply Harmony: ',
            lpos='w',
            wget_kwargs={'variable': self.harmony}
        ).grid(row=1, column=0, sticky='W')

        # - Embedding methods
        embedding_frame = ttk.Frame(self, borderwidth=1, relief='groove')
        embedding_frame.rowconfigure(0, weight=0)
        embedding_frame.rowconfigure(1, weight=1)
        embedding_frame.columnconfigure(0, weight=1)
        embedding_frame.columnconfigure(1, weight=1)
        embedding_frame.columnconfigure(2, weight=1)

        ttk.Label(
            embedding_frame,
            text='Select embedding method:'
        ).grid(row=0, columnspan=3, sticky='W')

        self.embedding_method = tk.StringVar()
        self.embedding_method.set('pca')
        ttk.Radiobutton(
            embedding_frame,
            text='PCA',
            variable=self.embedding_method,
            value='pca'
        ).grid(row=1, column=0, sticky='W')
        ttk.Radiobutton(
            embedding_frame,
            text='t-SNE',
            variable=self.embedding_method,
            value='tsne'
        ).grid(row=1, column=1, sticky='W')
        ttk.Radiobutton(
            embedding_frame,
            text='UMAP',
            variable=self.embedding_method,
            value='umap'
        ).grid(row=1, column=2, sticky='W')

        embedding_frame.grid(row=2, column=0, sticky='NSEW')

        # - Color variable selector
        self.color_var = tk.StringVar()
        self.combox_color_var = LabeledWidget(
            self,
            ttk.Combobox,
            'Select variable to color the embedding:',
            lpos='n',
            wget_kwargs={
                'textvariable': self.color_var,
                'state': 'disabled',
            },
            wget_grid_kwargs={'sticky': 'EW', 'weight': 1},
            label_grid_kwargs={'sticky': 'EW', 'weight': 0},
        )
        self.combox_color_var.grid(row=3, column=0, sticky='NSEW')
        self.combox_color_var.wg.bind('<<ComboboxSelected>>', self._update)

        ## PARAMS

        # - Plot button
        self.button_compute = ttk.Button(
            self,
            text='Compute',
            command=self.plot,
            state='disabled',
        )
        self.button_compute.grid(row=5, column=0, sticky='W')

        # Clustering panel
        ttk.Label(
            self,
            text='Clustering:',
            style='Title.TLabel',
            width=20
        ).grid(row=0, column=1, sticky='NSWE')

        # Plot
        self.fig, self.ax = plt.subplots()

        self.figframe = Figure(self, self.fig)
        self.figframe.grid(row=6, columnspan=2, sticky='NSWE')


    def _update(self, *ev):

        if self.controller.data and not self.controller.data.obs.empty:

            obs_key = self.color_var.get()

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
                    self.combox_color_var.wg.configure(
                        state='readonly',
                        values=obs_keys,
                    )
                    self.color_var.set(obs_key)

                self.button_compute.configure(state='normal')


    def plot(self):

        # TODO: apply Harmony?
        # TODO: plot_{embedding} should accept ax

        self.ax.clear()
        getattr(plots, f'plot_{self.embedding_method.get()}')(
            self.controller.data,
            color=self.color_var.get(),
            ax=self.ax,
        )
        self.figframe._update()