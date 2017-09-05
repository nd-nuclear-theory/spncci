#!/bin/csh

set archive_name = $1
set target_parent_dir = $2

set basename = ${archive_name:r}
set target_dir = ${target_parent_dir}/${basename}
echo "${archive_name} -> ${target_dir}"

echo "Target dir: ${target_dir}"
if (-e ${target_dir}) then
  echo "Target already exists."
else
  echo "Extracting..."
  mkdir --parents ${target_dir}
  nohup tar xf ${archive_name} --directory=${target_dir}
endif
