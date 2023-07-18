# ReaxPro wrappers:  

This repository is a collection of wrappers used in the [ReaxPro](https://www.reaxpro.eu/) project.
It holds wrappers for the [Amsterdam Modeling Suite](https://www.scm.com/amsterdam-modeling-suite/) (AMS) and [Zacros](https://zacros.org/).

## Authors
- [Pablo Lopez-Tarifa](mailto:p.lopez@esciencecenter.nl) (Main author). The Netherlands eScience Center.
- [Matthias Büschelberger](mailto:matthias.bueschelberger@iwm.fraunhofer.de) (Contributor). Fraunhofer Instituet for Mechanics of Materials (IWM)
- [Joana Francisco Morgado](mailto:joana.francisco.morgado@iwm.fraunhofer.de) (Contributor). Fraunhofer Institute for Mechanics of Materials (IWM)

## Index

- [Requirements](#requirements)
- [Wrapper generals](#wrapper-generals)
- [Installation](#installation)
- [Example](#example)

## Requiremnts

If you want to use this wrapper set, make sure you have installed:

- The [Simphony OSP-core](https://github.com/simphony/osp-core) version > 3.8.0.
- The [ReaxPro ontology](https://gitlab.cc-asp.fraunhofer.de/ontology/applications/reaxpro/reaxpro-framework).
- For AMS users, a licensed copy of [AMS](https://www.scm.com/amsterdam-modeling-suite/) installed.
- For Zacros users, both a licensed copy of [Zacros](https://zacros.org/software) code and [pyZacros](https://github.com/NLESC-JCER/pyZacros) library installed.

## Structure
A wrapper is a piece of code that slightly modifies the behavior of a function. 

The ReaxPro wrappers are built around the running functions of the above-mentioned software. Their main task is to translate (map) the semantic script provided by the user to the terms that are understood by the engines.

For a given engine XXX, there is a folder reaxpro-wrappers/osp/wrappers/simXXX/simXXX_session.py containing the Simphony wrapper session that will trigger the job execution.

In the folder reaxpro-wrappers/osp/tools is placed all the tooling for the semantic to syntactic mapping. 

## Installation
First of all, you will need to install OSP-core and pyZacros

```shell
(env) user@computer:~/reaxpro-wrappers$ pip install osp-core https://github.com/SCM-NV/pyZacros/archive/refs/tags/v.1.2.zip
```
Then, make sure that the wrapper can access the ontology from the [Fraunhofer Gitlab](https://gitlab.cc-asp.fraunhofer.de/) and download it with a given access token (with `read_api` and `read_repository` scopes). 

If you are using Windows, please type:

```shell
(env) C:\Users\user> set GITLAB_ACCESS_TOKEN=<your-access-token>
```

If you are using Linux, please type:

```shell
(env) user@computer:~/reaxpro-wrappers$ export GITLAB_ACCESS_TOKEN=<your-access-token>
```

Then, finally install the wrapper. Simply type:

```shell
(env) user@computer:~/reaxpro-wrappers$ python setup.py install
```

... or:


```shell
(env) user@computer:~/reaxpro-wrappers$ pip install .
```

## Example

The script [ams_wrapper.py](https://gitlab.cc-asp.fraunhofer.de/simphony/wrappers/reaxpro-wrappers/-/blob/master/examples/ams_wrapper.py) provides a simple semantic workflow to run a geometry optimization of a water molecule. 

To run the script:

```shell
(env) user@computer:~/reaxpro-wrappers$ python ams_wrapper.py 
```
