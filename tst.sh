#!/bin/bash

dstat -t -c --disk-util --top-cpu --top-io --top-mem -s --output system_logs.log 3 15 &

echo hi

