# add to end of path
append_path()
{
    # if zero length, just set equal to the input
    if eval test -z "\$$1"; then 
        eval "$1=$2"
    else
        if ! eval test -z "\"\${$1##*:$2:*}\"" -o -z "\"\${$1%%*:$2}\"" -o -z "\"\${$1##$2:*}\"" -o -z "\"\${$1##$2}\"" ; then
            eval "$1=\$$1:$2"
        fi
    fi
}

# add to front of path
prepend_path()
{
    # if zero length, just set equal to the input
    if eval test -z "\$$1"; then 
        eval "$1=$2"
    else
        if ! eval test -z "\"\${$1##*:$2:*}\"" -o -z "\"\${$1%%*:$2}\"" -o -z "\"\${$1##$2:*}\"" -o -z "\"\${$1##$2}\"" ; then
            eval "$1=$2:\$$1"
        fi
    fi
}



export PYORBFIT_HOME=/Users/gerdes/TNO/pyOrbfit

append_path LD_LIBRARY_PATH $PYORBFIT_HOME
append_path DYLD_LIBRARY_PATH $PYORBFIT_HOME
append_path C_INCLUDE_PATH $PYORBFIT_HOME
append_path PATH $PYORBFIT_HOME
append_path PYTHONPATH $PYORBFIT_HOME

export LD_LIBRARY_PATH 
export C_INCLUDE_PATH 
export PATH 

