#!/usr/bin/env nextflow

count = 0

process deseq2 {

    count++
    println count
    script:
        """
        echo 'Hello World'
        """
}
