# -*- coding: utf-8 -*-
"""
Created on Tue May 14 12:32:55 2019

@author: mmelkowski
"""


import time
import threading
import queue


class Thread_example:
    """
    The main idea behind this thread cnstruction is that you initialise a
    thread to perform a function. In this method the function need to be
    design to work once

    See Prime statistic for a working example.
    """

    def __init__(self, nb_thread_process=3):

        # NB of Thread for process
        self.nb_thread_process = nb_thread_process

        print("Start_main")
        self.main()
        print("Straight out of main")

    def worker(self):
        while True:
            time.sleep(1)
            print("in thread\n")
            if True:
                break

    def main(self):
        # Use of Lock to add to an object without having conflict
        # self.add_lock = threading.Lock()

        # Some object like queue.Queue will have integrated lock

        list_thread = []

        # THREAD init:
        for x in range(self.nb_thread_process):
            t = threading.Thread(target=self.worker) #args=() to pass argument in the function
            #  classifying as a daemon, so they will die when the main dies
            t.daemon = True
            #  begins, must come after daemon definition
            t.start()
            list_thread.append(t)

        # joining the thread will make the function wait until all thread are
        # finished before continuing
        for t in list_thread:
            t.join()


if __name__ == "__main__":
    start = time.time()
    stat = Thread_example(nb_thread_process=4)
    end = time.time() - start
    print("End in", round(end, 4), "sec")
