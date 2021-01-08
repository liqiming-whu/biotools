#!/usr/bin/env python3


def slice_by_len(length_list, container):
    results = []
    assert sum(length_list) == len(container), "length not match."
    for i, length in enumerate(length_list):
        start = sum(length_list[:i])
        end = start + length
        results.append(container[start:end])
    return results
