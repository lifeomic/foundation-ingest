#!groovy

pipeline {
  agent none

  stages {
    stage('build and security scan') {
      parallel {
        stage('build') {
          agent { label 'ecs-builder-node12' }
          steps {
            initBuild()
            sh 'yarn'
            sh 'yarn build'
          }
        }

        stage('security scan') {
          agent { label 'ecs-builder-node12' }
          steps {
            initBuild()
            sh 'yarn install'
            securityScan()
          }
        }
      }
    }
  }
}
