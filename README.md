### Xiayanwen Branch For Profiler


profiler to use :

### 1、gprof(both on master and slave) 

gprof工具的使用方法：在bsub脚本里加上-sw3runarg “-p -f”（主从核） -sw3runarg “-P master”（主核） -sw3runarg “-P slave”（从核）

http://bbs.nsccwx.cn/assets/uploads/files/1529047771770-%E7%A5%9E%E5%A8%81%E5%A4%AA%E6%B9%96%E4%B9%8B%E5%85%89%E9%AB%98%E7%BA%A7%E8%BF%9B%E9%98%B6.pdf

### 2、swlu(only for the master)

http://bbs.nsccwx.cn/topic/262/swlu-%E4%B8%BB%E6%A0%B8%E6%80%A7%E8%83%BD%E9%87%87%E6%A0%B7-%E8%B0%83%E8%AF%95%E5%B7%A5%E5%85%B7%E5%8C%85

### 3、lwpf(every kernel cache miss)

(can be used for slave)

http://bbs.nsccwx.cn/topic/257/lwpf2-%E8%BD%BB%E9%87%8F%E7%BA%A7%E4%BB%8E%E6%A0%B8%E6%80%A7%E8%83%BD%E6%8F%92%E6%A1%A9%E9%87%87%E6%A0%B7%E5%B7%A5%E5%85%B7

### 4、gptl (parallel)(seem not to be used for slave, or related docs cannot be found)

https://jmrosinski.github.io/GPTL/

