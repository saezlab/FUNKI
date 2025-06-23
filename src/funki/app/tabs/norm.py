import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt

from funki.pipelines import sc_quality_control
from funki.preprocessing import sc_trans_filter
from funki.preprocessing import sc_trans_normalize_total

from utils import Figure
from utils import LabeledWidget


class TabNorm(ttk.Frame):

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
        self.rowconfigure(5, weight=1)

        # Filter panel
        ttk.Label(
            self,
            text='Filtering:',
            style='Title.TLabel',
            width=20
        ).grid(row=0, column=0, sticky='NSWE')

        self.max_genes = tk.IntVar(value=5000)
        LabeledWidget(
            self,
            ttk.Entry,
            'Max. genes per cell: ',
            lpos='w',
            wget_kwargs={
                'textvariable': self.max_genes,
                'validate': 'key',
                'validatecommand': self.controller.check_num,
                'width': 10,
            }
        ).grid(row=1, column=0, sticky='W')

        self.min_genes = tk.IntVar(value=500)
        LabeledWidget(
            self,
            ttk.Entry,
            'Min. genes per cell: ',
            lpos='w',
            wget_kwargs={
                'textvariable': self.min_genes,
                'validate': 'key',
                'validatecommand': self.controller.check_num,
                'width': 10,
            }
        ).grid(row=2, column=0, sticky='W')

        self.mito_pct = tk.IntVar(value=5)
        LabeledWidget(
            self,
            tk.Scale,
            'Max. % mito. genes per cell: ',
            lpos='n',
            wget_kwargs={
                'from_': 0,
                'to': 100,
                'orient': 'horizontal',
                'variable': self.mito_pct,
            },
            wget_grid_kwargs={'sticky': 'WE'}
        ).grid(row=3, column=0, sticky='W')

        ttk.Button(
            self,
            text='Apply filters',
            command=self.apply_filter
        ).grid(row=4, column=0, sticky='W')

        # Normalization panel
        ttk.Label(
            self,
            text='Normalization:',
            style='Title.TLabel',
            width=20
        ).grid(row=0, column=1, sticky='NSWE')

        self.size_factor = tk.IntVar(value=1000000)
        LabeledWidget(
            self,
            ttk.Entry,
            'Size factor: ',
            lpos='w',
            wget_kwargs={
                'textvariable': self.size_factor,
                'validate': 'key',
                'validatecommand': self.controller.check_num,
                'width': 10,
            }
        ).grid(row=1, column=1, sticky='W')

        self.log_transform = tk.BooleanVar()
        LabeledWidget(
            self,
            ttk.Checkbutton,
            'Log-transform: ',
            lpos='w',
            wget_kwargs={'variable': self.log_transform}
        ).grid(row=2, column=1, sticky='W')

        ttk.Button(
            self,
            text='Apply normalization',
            command=self.apply_norm
        ).grid(row=4, column=1, sticky='W')

        # Figure
        self.fig, self.ax = plt.subplots(nrows=2, ncols=3)

        self.figframe = Figure(self, self.fig)
        self.figframe.grid(row=5, columnspan=2, sticky='NSWE')


    def _update(self, *ev):

        if self.controller.data:

            for ax in self.ax.flat:

                ax.clear()

            sc_quality_control(self.controller.data, ax=self.ax)
            self.figframe._update()


    def apply_filter(self, *ev):

        if self.controller.data:

            self.controller.data = sc_trans_filter(
                self.controller.data,
                min_genes=self.min_genes.get(),
                max_genes=self.max_genes.get(),
                mito_pct=self.mito_pct.get(),
            )

            self.controller._update()


    def apply_norm(self, *ev):

        if self.controller.data:

            self.controller.data = sc_trans_normalize_total(
                self.controller.data,
                target_sum=self.size_factor.get(),
                log_transform=self.log_transform.get()
            )

            self.controller._update()