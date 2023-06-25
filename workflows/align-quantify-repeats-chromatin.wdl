version 1.0

import "../tasks/star-smartmap.wdl" as smartmap

workflow wf {
    call smartmap.smartmap
}
