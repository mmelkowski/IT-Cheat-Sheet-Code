# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 12:17:16 2019

@author: mmelkowski
"""

# Exemple du fonctionnement d'asyncIO

import asyncio


async def count():
    print("One")
    await asyncio.sleep(1)
    print("Two")


async def main():
    await asyncio.gather(count(), count(), count())
