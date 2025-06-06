import tkinter as tk
from tkinter import ttk


def read_text(path):

    with open(path, 'r') as f:
    
        txt = ''.join(f.readlines())

    return txt


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

            # Wrap to the bounding box reserved for the label
            # instead of the container's full width.
            # Not doing this might lead to clipped text.
            options = self.grid_info()
            bbox = self.master.grid_bbox(options["column"], options["row"])

            if bbox is None:

                return 1

            width = bbox[2]

            if isinstance(options["padx"], int):

                padx = options["padx"] * 2

            else:

                padx = sum(options["padx"])

            return width - padx

        else:

            return self.master.winfo_width()