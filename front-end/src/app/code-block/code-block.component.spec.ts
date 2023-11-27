import { ComponentFixture, TestBed } from '@angular/core/testing';
import { By } from '@angular/platform-browser';
import { MatCardModule } from '@angular/material/card';
import { from } from 'rxjs';
import { CodeBlockComponent } from './code-block.component';
import { BlockService } from '../block.service';
import { HttpClientTestingModule } from '@angular/common/http/testing';

describe('CodeBlockComponent', () => {
  let component: CodeBlockComponent;
  let fixture: ComponentFixture<CodeBlockComponent>;

  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [CodeBlockComponent],
      imports: [
        HttpClientTestingModule,
        MatCardModule,
      ],
      providers: [BlockService],
    });
    fixture = TestBed.createComponent(CodeBlockComponent);
    component = fixture.componentInstance;
    component.block = {
      blockId: 'LoadData',
      title: 'Load Data',
      possibleChildBlocks: [],
      parameters: {},
      onRun: () => from(''),
    };
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });

  describe('RemoveBlock', () => {
    it ('should remove block when button is clicked', () => {
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('LoadData');
    });
  });
});
