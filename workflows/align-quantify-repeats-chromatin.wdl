version 1.0

import "../tasks/task-smartmap.wdl" as smartmap

workflow wf {
    call smartmap.smartmap
}
