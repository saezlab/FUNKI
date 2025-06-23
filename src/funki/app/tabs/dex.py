import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt

from funki.pipelines import differential_expression

from utils import Figure
from utils import LabeledWidget
from utils import PATH_ICON_SWP


class TabDex(ttk.Frame):

    def __init__ (self, parent, controller, **options):

        super().__init__(parent, **options)

        self.controller = controller

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=1)
        self.rowconfigure(0, weight=0)
        self.rowconfigure(1, weight=0)
        self.rowconfigure(2, weight=0)
        self.rowconfigure(3, weight=0)
        self.rowconfigure(4, weight=0)
        self.rowconfigure(5, weight=1)

        ttk.Label(
            self,
            text='Differential expression:',
            style='Title.TLabel',
            width=20
        ).grid(row=0, columnspan=3, sticky='NSWE')

        # Contrast variable selector
        self.contrast_var = tk.StringVar()
        self.combox_contrast_var = LabeledWidget(
            self,
            ttk.Combobox,
            'Select group variable for the contrast: ',
            lpos='n',
            wget_kwargs={
                'state': 'disabled',
                'textvariable': self.contrast_var,
            },
            wget_grid_kwargs={'sticky': 'EW', 'weight': 1},
            label_grid_kwargs={'sticky': 'EW', 'weight': 0},
        )
        self.combox_contrast_var.grid(row=1, columnspan=2, sticky='NSEW')
        self.combox_contrast_var.wg.bind('<<ComboboxSelected>>', self._update)

        # TODO: Handle multiple groups?

        # Contrast group A
        self.groupA = tk.StringVar()
        self.combox_A = LabeledWidget(
            self,
            ttk.Combobox,
            'Group of samples to contrast from: ',
            lpos='n',
            wget_kwargs={'state': 'disabled', 'textvariable': self.groupA},
            wget_grid_kwargs={'sticky': 'EW', 'weight': 1},
            label_grid_kwargs={'sticky': 'EW', 'weight': 0},
        )
        self.combox_A.grid(row=2, column=0, sticky='NSWE')
        self.combox_A.wg.bind('<<ComboboxSelected>>', self._update)

        # Contrast group B
        self.groupB = tk.StringVar()
        self.combox_B = LabeledWidget(
            self,
            ttk.Combobox,
            'Group of samples to contrast against: ',
            lpos='n',
            wget_kwargs={'state': 'disabled', 'textvariable': self.groupB},
            wget_grid_kwargs={'sticky': 'EW', 'weight': 1},
            label_grid_kwargs={'sticky': 'EW', 'weight': 0},
        )
        self.combox_B.grid(row=2, column=2, sticky='NSWE')
        self.combox_B.wg.bind('<<ComboboxSelected>>', self._update)

        # Swap button
        icon = tk.PhotoImage(file=PATH_ICON_SWP)
        self.button_swap = ttk.Button(
            self,
            image=icon,
            command=self.swap_ab,
            state='disabled',
        )
        self.button_swap.image = icon
        self.button_swap.grid(row=2, column=1, sticky='SW')

        # Thresholds
        self.thr_logfc = tk.DoubleVar()
        self.thr_logfc.set(1.0)
        LabeledWidget(
            self,
            ttk.Entry,
            'log(FC) threshold: ',
            lpos='w',
            wget_kwargs={
                'textvariable': self.thr_logfc,
                'validate': 'key',
                'validatecommand': self.controller.check_num,
                'width': 5,
            }
        ).grid(column=0, row=3, sticky='W')

        self.thr_pval = tk.DoubleVar()
        self.thr_pval.set(0.05)
        LabeledWidget(
            self,
            ttk.Entry,
            'P-value threshold: ',
            lpos='w',
            wget_kwargs={
                'textvariable': self.thr_pval,
                'validate': 'key',
                'validatecommand': self.controller.check_num,
                'width': 5,
            }
        ).grid(column=2, row=3, sticky='W')

        # DEX methods
        method_frame = ttk.Frame(
            self,
            borderwidth=1,
            relief='groove',
            padding=(5, 5, 5, 5),
        )
        method_frame.columnconfigure(0, weight=0)
        method_frame.columnconfigure(1, weight=1)
        method_frame.columnconfigure(2, weight=1)

        ttk.Label(
            method_frame,
            text='Select method:'
        ).grid(row=0, column=0, sticky='W')

        self.method = tk.StringVar()
        self.method.set('limma')
        ttk.Radiobutton(
            method_frame,
            text='limma',
            variable=self.method,
            value='limma'
        ).grid(row=0, column=1, sticky='E')
        ttk.Radiobutton(
            method_frame,
            text='PyDESeq2',
            variable=self.method,
            value='pydeseq2'
        ).grid(row=0, column=2, sticky='E')

        method_frame.grid(row=4, columnspan=2, sticky='NSEW')

        # Compute button
        self.button_compute = ttk.Button(
            self,
            text='Compute',
            command=self.compute,
            state='disabled',
        )
        self.button_compute.grid(row=4, column=2, sticky='SE')

        # Figure
        self.fig, self.ax = plt.subplots()

        self.figframe = Figure(self, self.fig)
        self.figframe.grid(row=5, columnspan=3, sticky='NSWE')


    def _update(self, *ev): # TODO: Could probably split to simplify

        if self.controller.data and not self.controller.data.obs.empty:

            obs_key = self.contrast_var.get()

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
                    self.contrast_var.set(obs_key)

            # Main combobox not empty
            if obs_key:

                options = sorted(set(self.controller.data.obs[obs_key]))

                # A/B are empty or new obs_key
                if (
                    (not self.groupA.get() or not self.groupB.wg.get())
                    or (
                        self.groupA.get() not in options
                        or self.groupB.get() not in options
                    )
                ):

                    self.combox_A.wg.configure(
                        state='readonly',
                        values=[n for i, n in enumerate(options) if i != 1],
                    )
                    self.groupA.set(options[0])

                    self.combox_B.wg.configure(
                        state='readonly',
                        values=[n for i, n in enumerate(options) if i != 0],
                    )
                    self.groupB.set(options[1])

                # A/B not empty -> ensure mutual exclusivity
                else:

                    curA = self.groupA.get()
                    curB = self.groupB.get()

                    self.combox_A.wg.configure(
                        values=[i for i in options if i != curB],
                    )

                    self.combox_B.wg.configure(
                        values=[i for i in options if i != curA],
                    )

                # Reactivating buttons
                self.button_compute.configure(state='normal')
                self.button_swap.configure(state='normal')


    def swap_ab(self):

        curA = self.groupA.get()
        curB = self.groupB.get()

        if curA and curB:

            obs_key = self.contrast_var.get()

            options = sorted(set(self.controller.data.obs[obs_key]))

            self.combox_A.wg.configure(
                values=[i for i in options if i != curA],
            )
            self.groupA.set(curB)
            self.combox_B.wg.configure(
                values=[i for i in options if i != curB],
            )
            self.groupB.set(curA)


    def compute(self):

        self.ax.clear()
        differential_expression(
            self.controller.data,
            self.contrast_var.get(),
            self.groupA.get(),
            self.groupB.get(),
            logfc_thr=self.thr_logfc.get(),
            fdr_thr=self.thr_pval.get(),
            method=self.method.get(),
            ax=self.ax,
        )
        self.figframe._update()