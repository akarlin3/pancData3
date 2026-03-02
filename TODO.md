---
trigger: always_on
---

1. Initialize the MATLAB environment by running 'source setup.sh'
2. Verify the license server variable is accessible: 'echo $MLM_LICENSE_FILE' (should return '27000@vsmphpllicmlab1')
3. Execute your code using the batch flag: 'matlab -batch "execute_all_workflows.m"'