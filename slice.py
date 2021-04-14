#!/usr/bin/env python3


def slice_by_len(length_list, container):
    container = list(container)
    container_len = len(container)
    length_sum = sum(length_list)
    results = []
    if container_len > length_sum:
        length_list.append(container_len - length_sum)
    for i, length in enumerate(length_list):
        start = sum(length_list[:i])
        end = start + length
        block = container[start:end]
        if block:
            results.append(block)
    return results
