#!/bin/bash
# Script based on http://tinyurl.com/lwvvgsn

set -o errexit -o nounset

if [ "$TRAVIS_BRANCH" != "master" ]
then
  echo "This commit was made against the $TRAVIS_BRANCH and not the master! No deploy!"
  exit 0
fi

if [ "$TRAVIS_PULL_REQUEST" == "true" ]
then
  echo "This pull request won't be deployed until it is accepted and merged! No deploy!"
  exit 0
fi

rev=$(git rev-parse --short HEAD)

cd doxygen/html

git init
git config user.name "Steven Gardiner"
git config user.email "sjgardiner@ucdavis.edu"

git remote add upstream "https://$GH_TOKEN@github.com/sjgardiner/marley.git"
git fetch upstream
git reset upstream/gh-pages

echo "www.marleygen.org" > CNAME

rm -f .gitignore
touch .

git add -A .
git commit -m "rebuild pages at ${rev}"
git push -q upstream HEAD:gh-pages
