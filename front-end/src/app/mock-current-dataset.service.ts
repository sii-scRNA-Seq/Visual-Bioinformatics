/* eslint-disable @typescript-eslint/no-unused-vars */

import { BehaviorSubject, Observable } from 'rxjs';
import { Injectable } from '@angular/core';

import { CurrentDatasetServiceInterface } from './current-dataset.service.interface';
import { DatasetInfo } from './dataset-info';

@Injectable({
  providedIn: 'root'
})
export class MockCurrentDatasetService implements CurrentDatasetServiceInterface {

  private readonly currentDataset$: BehaviorSubject<DatasetInfo> = new BehaviorSubject<DatasetInfo> ({
    key: '',
    title: '',
    samples: [],
    integration_obs: [],
    grouping_obs: []
  });
  readonly currentDataset: Observable<DatasetInfo> = this.currentDataset$.asObservable();

  setCurrentDataset(newDataset: string): void {
    this.currentDataset$.next({
      key: 'option1',
      title: 'Option 1',
      samples: ['sample1', 'sample2'],
      integration_obs: ['ob1', 'ob2'],
      grouping_obs: ['ob1', 'ob2']
    });
  }
}
