# ClonalTracker

Here is the ClonalTracker python code and databases

## Getting started

To make it easy for you to get started with GitLab, here's a list of recommended next steps.

Already a pro? Just edit this README.md and make it your own. Want to make it easy? [Use the template at the bottom](#editing-this-readme)!

## Add your files

- [ ] [Create](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#create-a-file) or [upload](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#upload-a-file) files
- [ ] [Add files using the command line](https://docs.gitlab.com/ee/gitlab-basics/add-file.html#add-a-file-using-the-command-line) or push an existing Git repository with the following command:

```
cd existing_repo
git remote add origin https://gitlab.com/victoriapascal/clonaltracker.git
git branch -M main
git push -uf origin main
```

## Integrate with your tools

- [ ] [Set up project integrations](https://gitlab.com/victoriapascal/clonaltracker/-/settings/integrations)

## Collaborate with your team

- [ ] [Invite team members and collaborators](https://docs.gitlab.com/ee/user/project/members/)
- [ ] [Create a new merge request](https://docs.gitlab.com/ee/user/project/merge_requests/creating_merge_requests.html)
- [ ] [Automatically close issues from merge requests](https://docs.gitlab.com/ee/user/project/issues/managing_issues.html#closing-issues-automatically)
- [ ] [Enable merge request approvals](https://docs.gitlab.com/ee/user/project/merge_requests/approvals/)
- [ ] [Automatically merge when pipeline succeeds](https://docs.gitlab.com/ee/user/project/merge_requests/merge_when_pipeline_succeeds.html)

## Test and Deploy

Use the built-in continuous integration in GitLab.

- [ ] [Get started with GitLab CI/CD](https://docs.gitlab.com/ee/ci/quick_start/index.html)
- [ ] [Analyze your code for known vulnerabilities with Static Application Security Testing(SAST)](https://docs.gitlab.com/ee/user/application_security/sast/)
- [ ] [Deploy to Kubernetes, Amazon EC2, or Amazon ECS using Auto Deploy](https://docs.gitlab.com/ee/topics/autodevops/requirements.html)
- [ ] [Use pull-based deployments for improved Kubernetes management](https://docs.gitlab.com/ee/user/clusters/agent/)
- [ ] [Set up protected environments](https://docs.gitlab.com/ee/ci/environments/protected_environments.html)

***

## How to run ClonalTracker
- Download ClonalTracker from gitlab (by using git clone or downloading a zip file): https://gitlab.com/victoriapascal/clonaltracker

```
git clone https://gitlab.com/victoriapascal/clonaltracker.git
```

- Install ClonalTracker dependencies by conda
```
conda env create -f clonaltracker/ct_env.yml ct_env
```

- Activate your environment
```
conda activate ct_env
```

- Run ClonalTracker
	- Fasta files should be in the same folder
	- Raw reads of both isolates should be in the same folder, and contain the *_R1.fastq.gz *_R2.fastq.gz suffix
```
python clonaltracker/clonaltracker.py /hpc/dla_mm/vpascalandreu/data/vanB_fastas/E7314.fasta /hpc/dla_mm/vpascalandreu/data/vanB_fastas/E7317.fasta /hpc/dla_mm/vpascalandreu/data/vanB_raw_reads_renamed2/
```

