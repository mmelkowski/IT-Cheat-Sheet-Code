# Multiprocessing

This folder describre the use of python module multiprocessing thorugh two main concept:

* pool
* queue

These template are taken from the following [Digital Ocean tutorial](https://www.digitalocean.com/community/tutorials/python-multiprocessing-example).

## Pool vs Queue

As [noxdafox](https://stackoverflow.com/a/48945440) says:

> The Pool and the Queue belong to two different levels of abstraction.
>
> The Pool of Workers is a concurrent design paradigm which aims to abstract a lot of logic you would otherwise need to implement yourself when using processes and queues.
>
> If your problem is simple enough, you can easily rely on a Pool. In more complex cases, you might need to deal with processes and queues yourself.
