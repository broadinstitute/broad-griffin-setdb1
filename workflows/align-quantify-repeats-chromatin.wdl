version 1.0

import "../tasks/task-smartmap.wdl" as smartmap

workflow wf {
    call smartmap.smartmap

    output{
        File smartmap_chromatin_prep = smartmap.smartmap_chromatin_prep
        File smartmap_bedgraph = smartmap.smartmap_bedgrap
    }

}
