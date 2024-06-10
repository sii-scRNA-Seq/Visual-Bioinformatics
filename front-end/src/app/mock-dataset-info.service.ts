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
      {key: 'option1', text: 'Option 1'},
      {key: 'option2', text: 'Option 2'}
    ]);
  }
}
