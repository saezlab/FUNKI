import re

import tkinter as tk
from tkinter import ttk
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk

from funki.common import _colors


PATH_LOGO = 'docs/source/_images/funki_logo.svg'


def read_text(path):

    with open(path, 'r') as f:
    
        txt = ''.join(f.readlines())

    return txt


def check_num(n):

    return re.match('^[0-9]*$', n) is not None


class LabeledWidget(ttk.Frame):
    '''
    Meta-widget consisting of a frame with a label in a specified position
    relative to a desired widget.
    '''

    def __init__(
        self,
        parent,
        Widget,
        label_text,
        lpos='N',
        wget_kwargs={},
        label_kwargs={},
        wget_grid_kwargs={},
        label_grid_kwargs={},
        **options
    ):

        super().__init__(parent, **options)

        self.wg = Widget(self, **wget_kwargs)
        self.label = ttk.Label(self, text=label_text, **label_kwargs)

        # Positioning
        lpos = lpos.lower()[0]

        if lpos not in 'nswe':

            raise ValueError(
                'Label position invalid, must be one of N, S, W or E.'
            )

        if lpos in 'ns':

            self.rowconfigure(0, weight=int(lpos == 's'))
            self.rowconfigure(1, weight=int(lpos == 'n'))

            self.wg.grid(column=0, row=int(lpos == 'n'), **wget_grid_kwargs)
            self.label.grid(column=0, row=int(lpos == 's'), **label_grid_kwargs)

        elif lpos in 'we':

            self.columnconfigure(0, weight=int(lpos == 'e'))
            self.columnconfigure(1, weight=int(lpos == 'w'))

            self.wg.grid(row=0, column=int(lpos == 'w'), **wget_grid_kwargs)
            self.label.grid(row=0, column=int(lpos == 'e'), **label_grid_kwargs)



class Figure(ttk.Frame):
    '''
    Container for a matplotlib figure with corresponding toolbar.
    '''

    def __init__(self, parent, fig, **options):

        super().__init__(parent, relief='sunken', borderwidth=1, **options)

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=0)
        
        self.canvas = FigureCanvasTkAgg(fig, master=self)

        toolbar = NavigationToolbar2Tk(self.canvas, self, pack_toolbar=False)
        toolbar.update()
        toolbar.grid(column=0, row=1, sticky='NSWE')
        
        self.canvas.get_tk_widget().grid(column=0, row=0, sticky='NSWE')

    
    def _update(self):

        self.canvas.draw()


class Table(ttk.Frame):
    '''
    Generates a ttk.Frame displaying a table with values extracted from a pandas
    DataFrame
    '''

    maxcols = 100
    maxrows = 100
    decimals = 3
    charlim = 25

    def __init__(self, parent, df=None, **options):

        super().__init__(parent, **options)

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=0)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=0)

        # Setting up
        self.tree = ttk.Treeview(self)
        self.tree.grid(column=0, row=0, sticky='NSEW')

        xscrollbar = ttk.Scrollbar(
            self,
            orient='horizontal',
            command=self.tree.xview,
        )
        xscrollbar.grid(column=0, row=1, sticky='EW')
        yscrollbar = ttk.Scrollbar(
            self,
            orient='vertical',
            command=self.tree.yview,
        )
        yscrollbar.grid(column=1, row=0, sticky='NS')

        self.tree.configure(
            xscrollcommand=xscrollbar.set,
            yscrollcommand=yscrollbar.set,
        )

        if df is not None:

            self.populate(df)

    
    def populate(self, df):
        '''
        Populates the table with the provided data.
        '''

        if not isinstance(df, pd.DataFrame):

            raise TypeError('Data provided is not a `pandas.DataFrame`')
        
        cropdf = self._crop_df(df)

        # Setting up contents
        columns = cropdf.columns.astype(str).to_list()
        self.tree['columns'] = columns
        self.tree.column('#0', width=150, anchor='w', stretch=False)

        for c in columns:

            self.tree.heading(c, text=c)
            self.tree.column(c, width=50, anchor='e', stretch=False)

        for i, (n, r) in enumerate(cropdf.iterrows()):

            self.tree.insert('', 'end', text=n, tag=i % 2, values=[
                self._fmt_val(v) for v in r
            ])

        self.tree.tag_configure(0, background=_colors['white'])
        self.tree.tag_configure(1, background=_colors['lightgray'])


    def _crop_df(self, df):
        '''
        Crops the data table based on maximum number of columns/rows for preview
        purposes.
        '''

        # Cropping rows
        if df.shape[0] > self.maxrows:

            top = bttm = self.maxrows // 2
            top += self.maxrows % 2

            x = pd.DataFrame(
                [['...'] * df.shape[1]],
                columns=df.columns,
                index=['...'],
            )

            df = pd.concat([df.iloc[:top, :], x, df.iloc[-bttm:, :]], axis=0)

        # Cropping columns
        if df.shape[1] > self.maxcols:

            left = right = self.maxcols // 2
            left += self.maxcols % 2

            x = pd.DataFrame(
                ['...'] * df.shape[0],
                columns=['...'],
                index=df.index,
            )

            df = pd.concat([df.iloc[:, :left], x, df.iloc[:, -right:]], axis=1)

        return df


    def _fmt_val(self, val):
        '''
        Formats a value in the table based on character/decimal limits. If the
        string is over the limit, cuts the characters to the limit and adds
        ellipsis. If the number has more decimals than the limit, it's rounded
        and cut up to those decimal positions.
        '''

        if isinstance(val, float):
            
            val = np.round(val, decimals=self.decimals)

        val = str(val)

        if len(val) > self.charlim:

            val = val[:self.charlim - 3] + '...'

        return val


# TODO: Should probably replace with Text widget, so we can add scrollbar and
#       also current solution does not wrap text very well
# Adapted from thegamecracks' gist
# https://gist.github.com/thegamecracks/5595dad631ec50bdc021f945054e86fb
class WrapLabel(ttk.Label):
    '''
    A label that automatically wraps its text when resized.

    When using this label, the geometry manager must be configured
    to resize the widget in some way, otherwise the width will
    default to the text's max length.
    This can mean using ``.pack(fill="both", expand=True)``, or
    ``.place(relwidth=1, relheight=1)``, or by configuring grid
    weights appropriately.

    ``minwidth=`` can be specified to prevent wrapping under
    a certain width which can significantly improve performance
    for longer lines of text.
    
    As a limitation of the current implementation, this class can only
    match its bounding box when used with the grid geometry manager.
    When using other geometry managers like pack or place, WrapLabel will
    always wrap to the window's entire width, regardless of the actual
    bounding box allocated for them. As a result, the text may appear clipped.
    
    If you use those geometry managers, it is recommended to grid / pack
    each WrapLabel inside a frame so it can calculate the wrap length from
    that frame's width. For example::
        
        # Instead of packing the label directly:
        label = WrapLabel(root, text="Some long piece of text")
        label.pack(fill="both", expand=True, padx=200, pady=200)
    
        # Pack a container frame and have the label expand inside it:
        frame = Frame(root)
        frame.pack(fill="both", expand=True, padx=200, pady=200)
        label = WrapLabel(frame, text="Some long piece of text")
        label.pack(fill="both", expand=True)
    
    In effect, the container frame serves as a way to compute the true
    bounding box allocated by the geometry manager, helping the label
    wrap to the correct width.
    '''

    def __init__(self, *args, minwidth: int = 1, **kwargs):

        super().__init__(*args, **kwargs)

        self.minwidth = minwidth
        self.bind("<Configure>", self.__on_configure)


    def __on_configure(self, event: tk.Event):

        width = max(self.minwidth, self.__get_width())

        if width != 1:  # Prevent wrapping on initial configuration

            self.configure(wraplength=width)


    def __get_width(self) -> int:

        if self.winfo_manager() == "grid":

            options = self.grid_info()
            col = int(options.get("column", 0))
            row = int(options.get("row", 0))
            colspan = int(options.get("columnspan", 1))

            # Get bbox for the full span (start_col to end_col)
            bbox = self.master.grid_bbox(col, row, col + colspan, row + 1)

            if bbox is None:

                return 1

            width = bbox[2]

            padx = options.get("padx", 0)

            if isinstance(padx, tuple):

                padx = sum(padx)

            else:

                padx *= 2

            return max(1, width - padx)

        else:

            return self.master.winfo_width()