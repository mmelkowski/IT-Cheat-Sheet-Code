# -*- coding: utf-8 -*-
"""
Created on Tue May 14 12:32:55 2019

@author: mmelkowski
"""

import tkinter as tk
# https://pythonprogramming.net/converting-tkinter-to-exe-with-cx-freeze/


class UI_example():
    def __init__(self):
        """
        In the positionning of the label we use .grid rather than .pack
        """
        # Creation of main frame
        self.root = tk.Tk()
        self.root.title("UI Main Title")

        #################################
        #            TITLE
        #################################
        i = 0

        label_Nom = tk.Label(self.root,
                             text="A Header",
                             font=("Consolas", 16))
        label_Nom.grid(row=i, column=1)

        #################################
        #           START GRID
        #################################
        i += 1

        # Label in place of the grid to show an introduction to the player
        intro = "\n  Text to show \n"
        self.label_center = tk.Label(self.root,
                                     justify=tk.LEFT,
                                     text=intro,
                                     font=("Consolas", 10))
        self.label_center.grid(row=i, column=1)

        # Label containing command legend for user
        legend = "Side text"
        label_side = tk.Label(self.root, text=legend)
        label_side.grid(row=i, column=2)

        #################################
        #        LINE FOR ENTRY
        #################################
        i += 1

        # Label in front of command input
        label_entree = tk.Label(self.root, text="Text to update:")
        label_entree.grid(row=i, column=0)

        # Entry widget to store command input
        self.entree = tk.Entry(self.root, textvariable="", width=30)
        self.entree.grid(row=i, column=1)

        # Button launch rover command execution
        button_Update = tk.Button(self.root,
                                  text=" Update Text ",
                                  command=self.update_text)
        button_Update.grid(row=i, column=2)

        #################################
        #           MAINLOOP
        #################################
        self.root.mainloop()

    def update_text(self):
        """
        Fonction of button: button_Update
        """
        # Get user inputs
        txt_in_label = self.entree.get()

        # Set New text for label_center
        self.label_center.configure(text=txt_in_label)
        # Refresh UI
        self.root.update_idletasks()


if __name__ == "__main__":
    UI_example()
