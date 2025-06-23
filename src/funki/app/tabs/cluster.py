import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt

from funki import plots
from funki.analysis import clustering
from funki.preprocessing import harmonize

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
        self.rowconfigure(6, weight=0)
        self.rowconfigure(7, weight=1)

        # Embedding panel
        ttk.Label(
            self,
            text='Embedding:',
            style='Title.TLabel',
            width=20
        ).grid(row=0, column=0, sticky='NSWE')

        # - Harmony
        self.harmony = tk.BooleanVar()
        LabeledWidget(
            self,
            ttk.Checkbutton,
            'Apply Harmony: ',
            lpos='w',
            wget_kwargs={
                'variable': self.harmony,
                'command': self._enable_harmony
            }
        ).grid(row=1, column=0, sticky='W')

        self.harmony_var = tk.StringVar()
        self.combox_harmony = LabeledWidget(
            self,
            ttk.Combobox,
            'Select variable to harmonize:',
            lpos='N',
            wget_kwargs={
                'state': 'disabled',
                'textvariable': self.harmony_var,
            },
            wget_grid_kwargs={'sticky': 'EW', 'weight': 1},
            label_grid_kwargs={'sticky': 'EW', 'weight': 0},
        )
        self.combox_harmony.grid(row=2, column=0, sticky='NSEW')
        self.combox_harmony.wg.bind('<<ComboboxSelected>>', self._update)

        # - Embedding methods
        embedding_frame = ttk.Frame(self)
        embedding_frame.rowconfigure(0, weight=0)
        embedding_frame.rowconfigure(1, weight=1)
        embedding_frame.columnconfigure(0, weight=1)
        embedding_frame.columnconfigure(1, weight=1)
        embedding_frame.columnconfigure(2, weight=1)

        ttk.Label(
            embedding_frame,
            text='Select embedding method:'
        ).grid(row=0, columnspan=3, sticky='W')

        self.embedding_method = tk.StringVar(value='pca')
        self.embedding_method_buttons = {
            v: ttk.Radiobutton(
                embedding_frame,
                text=t,
                variable=self.embedding_method,
                command=self.set_embed_param,
                value=v,
            )
            for t, v in [('PCA', 'pca'), ('t-SNE', 'tsne'), ('UMAP', 'umap')]
        }

        for i, wg in enumerate(self.embedding_method_buttons.values()):

            wg.grid(row=1, column=i, sticky='W')

        embedding_frame.grid(row=3, column=0, sticky='NSEW')

        # - Embedding parameters frame (contents set by set_embed_param)
        self.perplexity = tk.DoubleVar(value=30)
        self.min_dist = tk.DoubleVar(value=0.5)
        self.spread = tk.DoubleVar(value=1.0)
        self.alpha = tk.DoubleVar(value=1.0)
        self.gamma = tk.DoubleVar(value=1.0)
        self.embedding_params_frame = ttk.Frame(
            self,
            borderwidth=1,
            relief='groove',
            padding=(5, 5, 5, 5),
        )
        self.embedding_params_frame.columnconfigure(0, weight=1)
        self.embedding_params_frame.columnconfigure(1, weight=1)
        self.embedding_params_frame.rowconfigure(0, weight=0)
        self.embedding_params_frame.rowconfigure(1, weight=1)
        self.embedding_params_frame.rowconfigure(2, weight=1)

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
        self.combox_color_var.grid(row=5, column=0, sticky='NSEW')
        self.combox_color_var.wg.bind('<<ComboboxSelected>>', self._update)

        # - Plot button
        self.button_plot = ttk.Button(
            self,
            text='Plot embedding',
            command=self.plot,
            state='disabled',
        )
        self.button_plot.grid(row=6, column=0, sticky='W')

        # Clustering panel
        ttk.Label(
            self,
            text='Clustering:',
            style='Title.TLabel',
            width=20
        ).grid(row=0, column=1, sticky='NSWE')

        # - Resoulution
        self.resoultion = tk.DoubleVar(value=1.0)
        LabeledWidget(
            self,
            ttk.Entry,
            'Resolution: ',
            lpos='w',
            wget_kwargs={
                'textvariable': self.resoultion,
                'validate': 'key',
                'validatecommand': self.controller.check_num,
                'width': 5,
            }
        ).grid(row=1, column=1, sticky='NSWE')

        # - Choice algorithm
        clustering_frame = ttk.Frame(self)
        clustering_frame.columnconfigure(0, weight=1)
        clustering_frame.columnconfigure(1, weight=1)
        clustering_frame.rowconfigure(0, weight=0)
        clustering_frame.rowconfigure(1, weight=1)
        clustering_frame.grid(row=2, column=1, sticky='NSWE')

        ttk.Label(
            clustering_frame,
            text='Select clustering method:'
        ).grid(row=0, columnspan=2, sticky='W')

        self.clustering_method = tk.StringVar(value='leiden')
        ttk.Radiobutton(
                clustering_frame,
                text='Leiden',
                variable=self.clustering_method,
                value='leiden'
        ).grid(row=1, column=0, sticky='W')
        ttk.Radiobutton(
                clustering_frame,
                text='Louvain',
                variable=self.clustering_method,
                value='louvain'
        ).grid(row=1, column=1, sticky='W')

        # - Compute button
        self.button_compute = ttk.Button(
            self,
            text='Compute',
            command=self.cluster,
            state='disabled',
        )
        self.button_compute.grid(row=6, column=1, sticky='W')

        # Plot
        self.fig, self.ax = plt.subplots()

        self.figframe = Figure(self, self.fig)
        self.figframe.grid(row=7, columnspan=2, sticky='NSWE')


    def _update(self, *ev):

        obs_key_color = self.color_var.get()
        obs_key_harmony = self.harmony_var.get()

        if self.controller.data and not self.controller.data.obs.empty:

            obs_keys = sorted(self.controller.data.obs_keys())

            if obs_keys:

                obs_key_color = obs_key_color or obs_keys[0]
                obs_key_harmony = obs_key_harmony or obs_keys[0]

                self.combox_color_var.wg.configure(
                    state='readonly',
                    values=obs_keys,
                )
                self.combox_harmony.wg.configure(
                    state='readonly' if self.harmony.get() else 'disabled',
                    values=obs_keys,
                )

                self.color_var.set(obs_key_color)
                self.harmony_var.set(obs_key_harmony)

            self.button_plot.configure(state='normal')
            self.button_compute.configure(state='normal')


    def _enable_harmony(self, *ev):

        state = (
            'readonly'
            if (self.harmony.get() and self.harmony_var.get())
            else 'disabled'
        )

        self.combox_harmony.wg.configure(state=state)


    def set_embed_param(self, *ev):

        # Clearing frame
        for child in self.embedding_params_frame.winfo_children():

            child.grid_forget()

        self.embedding_params_frame.grid_forget()

        method = self.embedding_method.get()

        if method != 'pca':

            self.embedding_params_frame.grid(
                row=4,
                column=0,
                sticky='NWE',
                padx=(10, 10),
            )

            mlabel = self.embedding_method_buttons[method].cget('text')
            ttk.Label(
                self.embedding_params_frame,
                text=f'{mlabel} parameters:',
            ).grid(row=0, columnspan=2, sticky='W')


        if method == 'tsne':

            LabeledWidget(
                self.embedding_params_frame,
                ttk.Entry,
                'Perplexity: ',
                lpos='w',
                wget_kwargs={
                    'textvariable': self.perplexity,
                    'validate': 'key',
                    'validatecommand': self.controller.check_num,
                    'width': 5,
                }
            ).grid(row=1, columnspan=2, sticky='W')

        elif method == 'umap':

            LabeledWidget(
                self.embedding_params_frame,
                ttk.Entry,
                'Min. distance: ',
                lpos='w',
                wget_kwargs={
                    'textvariable': self.min_dist,
                    'validate': 'key',
                    'validatecommand': self.controller.check_num,
                    'width': 5,
                }
            ).grid(row=1, column=0, sticky='W')
            LabeledWidget(
                self.embedding_params_frame,
                ttk.Entry,
                'Spread: ',
                lpos='w',
                wget_kwargs={
                    'textvariable': self.spread,
                    'validate': 'key',
                    'validatecommand': self.controller.check_num,
                    'width': 5,
                }
            ).grid(row=1, column=1, sticky='W')
            LabeledWidget(
                self.embedding_params_frame,
                ttk.Entry,
                'Alpha: ',
                lpos='w',
                wget_kwargs={
                    'textvariable': self.alpha,
                    'validate': 'key',
                    'validatecommand': self.controller.check_num,
                    'width': 5,
                }
            ).grid(row=2, column=0, sticky='W')
            LabeledWidget(
                self.embedding_params_frame,
                ttk.Entry,
                'Gamma: ',
                lpos='w',
                wget_kwargs={
                    'textvariable': self.gamma,
                    'validate': 'key',
                    'validatecommand': self.controller.check_num,
                    'width': 5,
                }
            ).grid(row=2, column=1, sticky='W')


    def cluster(self, *ev):

        clustering(
            self.controller.data,
            alg=self.clustering_method.get(),
            resolution=self.resoultion.get(),
        )
        self._update()


    def plot(self):

        # Apply Harmony?
        if self.harmony.get() and self.harmony_var.get():

            harmonize(
                self.controller.data,
                [self.harmony_var.get()],
                recalculate=True
            )

        method = self.embedding_method.get()

        kwargs = {}

        if method == 'tsne':

            kwargs['perplexity'] = self.perplexity.get()

        elif method == 'umap':

            kwargs['min_dist'] = self.min_dist.get()
            kwargs['spread'] = self.spread.get()
            kwargs['alpha'] = self.alpha.get()
            kwargs['gamma'] = self.gamma.get()

        self.ax.clear()
        getattr(plots, f'plot_{method}')(
            self.controller.data,
            color=self.color_var.get(),
            ax=self.ax,
            **kwargs
        )
        self.figframe._update()
