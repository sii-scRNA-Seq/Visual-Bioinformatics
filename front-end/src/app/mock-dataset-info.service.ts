import { BehaviorSubject, Observable } from 'rxjs';
import { Injectable } from '@angular/core';

import { DatasetInfo } from './dataset-info';
import { DatasetInfoServiceInterface } from './dataset-info.service.interface';

@Injectable({
  providedIn: 'root'
})
export class MockDatasetInfoService implements DatasetInfoServiceInterface {

  private readonly datasetInfo$: BehaviorSubject<DatasetInfo[]> = new BehaviorSubject<DatasetInfo[]> ([]);
  readonly datasetInfo: Observable<DatasetInfo[]> = this.datasetInfo$.asObservable();

  async setDatasetInfo(): Promise<void> {
    this.datasetInfo$.next([
      {
        key: 'option1',
        title: 'Option 1',
        samples: ['sample1', 'sample2'],
        integration_obs: ['ob1', 'ob2'],
        grouping_obs: ['ob1', 'ob2']
      },
      {
        key: 'option2',
        title: 'Option 2',
        samples: [],
        integration_obs: [],
        grouping_obs: []
      }
    ]);
  }
}
